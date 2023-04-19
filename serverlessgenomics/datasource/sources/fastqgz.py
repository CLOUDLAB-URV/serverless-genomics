from __future__ import annotations

import io
import logging
import os
import re
import subprocess
import tempfile
import threading
import time
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from serverlessgenomics.pipeline import PipelineParameters, Lithops
from serverlessgenomics.utils import try_head_object, S3Path, force_delete_local_path

if TYPE_CHECKING:
    from typing import Tuple, List
    from lithops import Storage

logger = logging.getLogger(__name__)

CHUNK_SIZE = 65536
RE_WINDOWS = re.compile(r"#\d+: @ \d+ / \d+ L\d+ \( \d+ @\d+ \)")
RE_NUMS = re.compile(r"\d+")
RE_NLINES = re.compile(r"Number of lines\s+:\s+\d+")


def check_fastqgz_index(pipeline_params: PipelineParameters, lithops: Lithops) -> int:
    """
    Check if gzip index file for FASTQ in storage already exists, or create one if it doesn't.
    Returns number of lines of FASTQ file.
    """
    # Check if fastqgz file exists
    fastq_head = try_head_object(lithops.storage, pipeline_params.fastq_path.bucket, pipeline_params.fastq_path.key)
    if fastq_head is None:
        raise Exception(f"FASTQGZip file with key {pipeline_params.fastq_path} does not exists")

    # Check if fastqgz index file exists
    index_key, tab_key = get_fastqgz_idx_keys(pipeline_params)
    fastq_idx_head = try_head_object(lithops.storage, pipeline_params.fastq_path.bucket, index_key)
    fastq_tab_head = try_head_object(lithops.storage, pipeline_params.fastq_path.bucket, tab_key)
    if None in (fastq_idx_head, fastq_tab_head):
        # Generate gzip index file for compressed fastq input
        logger.info("Generating gzip index file for FASTQ %s", pipeline_params.fastq_path.stem)
        total_lines, dictio = lithops.invoker.call(
            generate_idx_from_gzip, (pipeline_params, pipeline_params.fastq_path)
        )
    else:
        # Get total lines from header metadata
        logger.debug("FASTQGZip index for %s found", pipeline_params.fastq_path.stem)
        total_lines = int(fastq_tab_head["x-amz-meta-total_lines"])

    return total_lines


def generate_idx_from_gzip(pipeline_params: PipelineParameters, gzip_file_path: S3Path, storage: Storage):
    """
    Lithops callee function
    Create index file from gzip archive using gztool (https://github.com/circulosmeos/gztool)
    """
    gztool = get_gztool_path()
    tmp_index_file_name = tempfile.mktemp()
    gzip_idx_key, gzip_tab_key = get_fastqgz_idx_keys(pipeline_params)

    try:
        data_stream = storage.get_object(bucket=gzip_file_path.bucket, key=gzip_file_path.key, stream=True)

        force_delete_local_path(tmp_index_file_name)
        t0 = time.perf_counter()

        # Create index and save to tmp file
        # TODO tmp file is needed, sending to stdout is not working at the moment (todo fix)
        index_proc = subprocess.Popen(
            [gztool, "-i", "-x", "-I", tmp_index_file_name],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        # TODO program might get stuck if subprocess fails, blocking io should be done in a backgroun thread or using
        #  async/await
        try:
            chunk = data_stream.read(CHUNK_SIZE)
            while chunk != b"":
                # logger.debug('Writing %d bytes to Pipe STDIN', len(chunk))
                index_proc.stdin.write(chunk)
                chunk = data_stream.read(CHUNK_SIZE)
            if hasattr(data_stream, "close"):
                data_stream.close()
        except BrokenPipeError as e:
            stdout, stderr = index_proc.communicate()
            logger.error(stdout.decode("utf-8"))
            logger.error(stderr.decode("utf-8"))
            raise e

        stdout, stderr = index_proc.communicate()
        logger.debug(stdout.decode("utf-8"))
        logger.debug(stderr.decode("utf-8"))
        if index_proc.returncode > 0:
            logger.debug(stdout.decode("utf-8"))
            logger.debug(stderr.decode("utf-8"))
            raise Exception("Error creating gz index")

        # Generate list of access windows from index
        proc = subprocess.run(
            [gztool, "-ell", "-I", tmp_index_file_name],
            check=True,
            capture_output=True,
            text=True,
        )
        output = proc.stdout
        # logger.debug(output)

        # Store index binary file
        storage.upload_file(
            file_name=tmp_index_file_name,
            bucket=pipeline_params.storage_bucket,
            key=gzip_idx_key,
        )

        # Get the total number of lines
        total_lines = int(RE_NUMS.findall(RE_NLINES.findall(output).pop()).pop())
        logger.debug("Indexed gzipped text file with %s total lines", total_lines)
        t1 = time.perf_counter()
        logger.debug("Index generated in %.3f seconds", t1 - t0)

        # Generator function that parses output to avoid copying all window data as lists
        def _lines_generator():
            for f in RE_WINDOWS.finditer(output):
                nums = [int(n) for n in RE_NUMS.findall(f.group())]
                yield nums

        # Generate data frame that stores gzip index windows offsets
        df = pd.DataFrame(
            _lines_generator(),
            columns=[
                "window",
                "compressed_byte",
                "uncompressed_byte",
                "line_number",
                "window_size",
                "window_offset",
            ],
        )
        df.set_index(["window"], inplace=True)

        # Store data frame as parquet
        out_stream = io.BytesIO()
        df.to_parquet(out_stream, engine="pyarrow")
        # df.to_csv(os.path.join(tempfile.gettempdir(), f'{meta.obj_path.stem}.csv'))  # debug

        os.remove(tmp_index_file_name)

        out_stream.seek(0)
        storage.get_client().upload_fileobj(
            Bucket=pipeline_params.storage_bucket,
            Key=gzip_tab_key,
            Fileobj=out_stream,
            ExtraArgs={"Metadata": {"total_lines": str(total_lines)}},
        )
    finally:
        force_delete_local_path(tmp_index_file_name)


def get_ranges_from_line_pairs(pipeline_params: PipelineParameters, pairs: List[Tuple[int, int]], storage: Storage):
    _, gzip_tab_key = get_fastqgz_idx_keys(pipeline_params)

    # Download gzip tab into an in-memory buffer and read dataframe into Pandas
    buff = io.BytesIO()
    storage.get_client().download_fileobj(Bucket=pipeline_params.storage_bucket, Key=gzip_tab_key, Fileobj=buff)
    buff.seek(0)
    df = pd.read_parquet(buff)
    del buff

    # Get fastq file size, required to calculate maximum offset for last chunk
    head_fastq = storage.head_object(bucket=pipeline_params.fastq_path.bucket, key=pipeline_params.fastq_path.key)
    fastqgz_sz = int(head_fastq["content-length"])

    line_indexes = df["line_number"].to_numpy()
    num_windows = df.shape[0]

    byte_ranges = [None] * len(pairs)
    for i, (line_0, line_1) in enumerate(pairs):
        # Find the closest window index for line_0
        window_head_idx = (np.abs(line_indexes - line_0)).argmin()
        # Check if window line entry pont is past requested line_0, if so, get previous window
        window_head_line = df.iloc[window_head_idx]["line_number"]
        if window_head_line > line_0:
            window_head_idx = window_head_idx - 1
        # Get offset in compressed archive for window 0
        window0_offset = df.iloc[window_head_idx]["compressed_byte"]

        # Find the closest window index for line_0
        widow_tail_idx = (np.abs(line_indexes - line_1)).argmin()
        # Check if window line entry pont is before requested line_1, if so, get next window
        window_tail_line = df.iloc[widow_tail_idx]["line_number"]
        if window_tail_line < line_1:
            widow_tail_idx = widow_tail_idx + 1
        if widow_tail_idx >= num_windows:
            # Adjust offset for lines inside last window, use end of compressed archive for 2nd offset
            window1_offset = fastqgz_sz
        else:
            window1_offset = df.iloc[widow_tail_idx]["compressed_byte"]

        byte_ranges[i] = (window0_offset, window1_offset)

    return byte_ranges


def fetch_fastq_chunk_s3_fastqgzip(
    fastq_chunk: dict, target_filename: str, pipeline_parameters: PipelineParameters, storage: Storage
):
    tmp_index_file = tempfile.mktemp()
    _, gzip_idx_key = get_fastqgz_idx_keys(pipeline_parameters)
    gztool = get_gztool_path()
    lines = []
    lines_to_read = fastq_chunk["line_1"] - fastq_chunk["line_0"] + 1

    try:
        t0 = time.perf_counter()

        # Get index and store it to temp file
        storage.download_file(bucket=pipeline_parameters.storage_bucket, key=gzip_idx_key, file_name=tmp_index_file)

        # Get compressed byte range
        extra_get_args = {"Range": f"bytes={fastq_chunk['range_0'] - 1}-{fastq_chunk['range_1'] - 1}"}
        body = storage.get_object(
            pipeline_parameters.fastq_path.bucket, pipeline_parameters.fastq_path.key, True, extra_get_args
        )

        cmd = [gztool, "-I", tmp_index_file, "-n", str(fastq_chunk["range_0"]), "-L", str(fastq_chunk["line_0"])]
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

        # TODO program might get stuck if subprocess fails, blocking io should be done in a backgroun thread or using async/await
        def _writer_feeder():
            logger.debug("Writer thread started")
            input_chunk = body.read(CHUNK_SIZE)
            while input_chunk != b"":
                # logger.debug('Writing %d bytes to pipe', len(chunk))
                try:
                    proc.stdin.write(input_chunk)
                except BrokenPipeError:
                    break
                input_chunk = body.read(CHUNK_SIZE)
            try:
                proc.stdin.flush()
                proc.stdin.close()
            except BrokenPipeError:
                pass
            logger.debug("Writer thread finished")

        writer_thread = threading.Thread(target=_writer_feeder)
        writer_thread.start()

        output_chunk = proc.stdout.read(CHUNK_SIZE)
        last_line = None
        while output_chunk != b"":
            # logger.debug('Read %d bytes from pipe', len(chunk))
            text = output_chunk.decode("utf-8")
            chunk_lines = text.splitlines()
            if last_line is not None:
                last_line = last_line + chunk_lines.pop(0)
                lines.append(last_line)
                last_line = None
            if text[-1] != "\n":
                last_line = chunk_lines.pop()

            lines.extend(chunk_lines)

            # Stop decompressing lines if number of lines to read in this chunk is reached
            if len(lines) > lines_to_read:
                proc.stdout.close()
                break

            # Try to read next decompressed chunk
            # a ValueError is raised if the pipe is closed, meaning the writer or the subprocess closed it
            try:
                output_chunk = proc.stdout.read(CHUNK_SIZE)
            except ValueError:
                output_chunk = b""

        try:
            proc.wait()
        except ValueError as e:
            logger.error(e)

        writer_thread.join()

        t1 = time.perf_counter()
        logger.debug("Got partition in %.3f seconds", t1 - t0)
        # TODO write lines to file as decompressed instead of saving them all in memory
        with open(target_filename, "w") as target_file:
            target_file.writelines((line + "\n" for line in lines[: fastq_chunk["line_1"] - fastq_chunk["line_0"]]))
    finally:
        force_delete_local_path(tmp_index_file)


def get_fastqgz_idx_keys(pipeline_params: PipelineParameters) -> Tuple[str, str]:
    """
    Helper function to format fastqgz index keys in storage as tuple (index file key, tab data key)
    """
    return (
        os.path.join(pipeline_params.fastqgz_idx_prefix, pipeline_params.fastq_path.key + ".idx"),
        os.path.join(pipeline_params.fastqgz_idx_prefix, pipeline_params.fastq_path.key + ".tab"),
    )


def get_gztool_path():
    """
    Utility function that returns the absolute path for gzip file binary or raises exception if it is not found
    """
    proc = subprocess.run(["which", "gztool"], check=True, capture_output=True, text=True)
    path = proc.stdout.rstrip("\n")
    logger.debug("Using gztool located in %s", path)
    return path


#
