from __future__ import annotations

from ..stats import Stats

import io
import logging
import re
import subprocess
import tempfile
import time
from math import ceil
from typing import TYPE_CHECKING, List, Tuple
import requests
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import requests

from ..utils import force_delete_local_path, S3Path, try_head_object, get_gztool_path

if TYPE_CHECKING:
    from ..parameters import PipelineRun, Lithops
    import lithops

logger = logging.getLogger(__name__)

CHUNK_SIZE = 65536
RE_WINDOWS = re.compile(r"#\d+: @ \d+ / \d+ L\d+ \( \d+ @\d+ \)")
RE_NUMS = re.compile(r"\d+")
RE_NLINES = re.compile(r"Number of lines\s+:\s+\d+")


# @lithops_callee
def generate_idx_from_gzip(pipeline_params: PipelineRun, gzip_file_path: S3Path, storage: lithops.Storage):
    """
    Lithops callee function
    Create index file from gzip archive using gztool (https://github.com/circulosmeos/gztool)
    """
    gztool = get_gztool_path()
    tmp_index_file_name = tempfile.mktemp()
    s3 = storage.get_client()

    stats = Stats()
    try:
        stats.timer_start("get_fastq_file")
        res = s3.get_object(Bucket=gzip_file_path.bucket, Key=gzip_file_path.key)
        stats.timer_stop("get_fastq_file")
        stats.store_size_data(
            "size", storage.head_object(gzip_file_path.bucket, gzip_file_path.key)["content-length"], "fastq_file"
        )

        data_stream = res["Body"]
        force_delete_local_path(tmp_index_file_name)
        t0 = time.perf_counter()

        # Create index and save to tmp file
        index_proc = subprocess.Popen(
            [gztool, "-i", "-x", "-I", tmp_index_file_name],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        chunk = data_stream.read(CHUNK_SIZE)
        while chunk != b"":
            index_proc.stdin.write(chunk)
            chunk = data_stream.read(CHUNK_SIZE)
        if hasattr(data_stream, "close"):
            data_stream.close()

        stdout, stderr = index_proc.communicate()
        # logger.debug(stdout.decode('utf-8'))
        # logger.debug(stderr.decode('utf-8'))
        if index_proc.returncode > 0:
            logger.debug(stdout.decode("utf-8"))
            logger.debug(stderr.decode("utf-8"))
            raise Exception("Error creating gz index")

        # Generate list of access windows from index
        proc = subprocess.run([gztool, "-ell", "-I", tmp_index_file_name], check=True, capture_output=True, text=True)
        output = proc.stdout

        # Store index binary file
        index_key, tab_key = pipeline_params.fastqgz_idx_keys
        stats.timer_start("upload_index_binary_fastq")
        s3.upload_file(Filename=tmp_index_file_name, Bucket=pipeline_params.storage_bucket, Key=index_key)
        stats.timer_stop("upload_index_binary_fastq")
        stats.store_size_data(
            "size",
            storage.head_object(bucket=pipeline_params.storage_bucket, key=index_key)["content-length"],
            "index_binary_fastq",
        )

        # Get the total number of lines
        total_lines: int = int(RE_NUMS.findall(RE_NLINES.findall(output).pop()).pop())
        logger.debug("Indexed gzipped text file with %d total lines", total_lines)
        t1 = time.perf_counter()
        logger.debug("Index generated in %.3f seconds", t1 - t0)

        # Generator function that parses output to avoid copying all window data as lists
        def _lines_generator():
            for f in RE_WINDOWS.finditer(output):
                # nums = np.array([int(n) for n in RE_NUMS.findall(f.group())], dtype=np.int)
                nums = [int(n) for n in RE_NUMS.findall(f.group())]
                yield nums

        # Generate data frame that stores gzip index windows offsets

        df = pd.DataFrame(
            _lines_generator(),
            columns=["window", "compressed_byte", "uncompressed_byte", "line_number", "window_size", "window_offset"],
        )
        df.set_index(["window"], inplace=True)

        # Store data frame as parquet
        out_stream = io.BytesIO()
        df.to_parquet(out_stream, engine="pyarrow")
        out_stream.seek(0)
        stats.timer_start("put_data_frame_parquet")
        s3.put_object(
            Bucket=pipeline_params.storage_bucket,
            Key=tab_key,
            Body=out_stream,
            Metadata={"total_lines": str(total_lines)},
        )
        stats.timer_stop("put_data_frame_parquet")
        stats.store_size_data(
            "size",
            storage.head_object(pipeline_params.storage_bucket, tab_key)["content-length"],
            "put_data_frame_parquet",
        )

        return total_lines, stats.get_stats()
    finally:
        force_delete_local_path(tmp_index_file_name)


def get_ranges_from_line_pairs(pipeline_params: PipelineRun, lithops: Lithops, pairs: List[Tuple[int, int]], stats):
    _, tab_key = pipeline_params.fastqgz_idx_keys
    stats.timer_start("get_data_frame_parquet")
    meta_obj_body = lithops.storage.get_object(bucket=pipeline_params.storage_bucket, key=tab_key)
    stats.timer_stop("get_data_frame_parquet")
    stats.store_size_data(
        "size",
        lithops.storage.head_object(pipeline_params.storage_bucket, tab_key)["content-length"],
        "get_data_frame_parquet",
    )

    meta_buff = io.BytesIO(meta_obj_body)
    meta_buff.seek(0)
    df = pd.read_parquet(meta_buff)
    line_indexes = df["line_number"].to_numpy()
    num_windows = df.shape[0]

    head_fastq = lithops.storage.head_object(
        bucket=pipeline_params.fastq_path.bucket, key=pipeline_params.fastq_path.key
    )
    fastqgz_sz = int(head_fastq["content-length"])

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


def generate_fastqgz_index_from_s3(pipeline_params: PipelineRun, lithops: Lithops) -> int:
    """
    Generate gzip index file for FASTQ in storage if necessary, returns number of lines of FASTQ file
    """
    # Check if fastq file exists
    fastq_head = try_head_object(lithops.storage, pipeline_params.fastq_path.bucket, pipeline_params.fastq_path.key)
    if fastq_head is None:
        raise Exception(f"fastq file with key {pipeline_params.fastq_path} does not exists")

    # Check if fastqgz index file exists
    index_key, tab_key = pipeline_params.fastqgz_idx_keys
    fastq_idx_head = try_head_object(lithops.storage, pipeline_params.fastq_path.bucket, index_key)
    fastq_tab_head = try_head_object(lithops.storage, pipeline_params.fastq_path.bucket, tab_key)
    dictio = {}
    if fastq_idx_head is None or fastq_tab_head is None:
        # Generate gzip index file for compressed fastq input
        logger.info("Generating gzip index file for FASTQ %s", pipeline_params.fastq_path.stem)
        res = lithops.invoker.call(generate_idx_from_gzip, (pipeline_params, pipeline_params.fastq_path))
        total_lines = res[0]
        dictio = res[1]
    else:
        # Get total lines from header metadata
        logger.debug("Fastqz index for %s found", pipeline_params.fastq_path.stem)
        total_lines = int(fastq_tab_head["x-amz-meta-total_lines"])
    logger.info("Read %d sequences from FASTQ %s", total_lines / 4, pipeline_params.fastq_path.stem)

    return total_lines, dictio


def prepare_fastq_chunks(pipeline_params: PipelineRun, lithops: Lithops):
    """
    Calculate fastq byte ranges for chunks of a pipeline run, generate gzip index if needed
    """
    subStat = Stats()
    subStat.timer_start("prepare_fastq_chunks")
    if pipeline_params.fastq_path is not None:
        res = generate_fastqgz_index_from_s3(pipeline_params, lithops)
        num_lines = res[0]
        subStat.store_dictio(res[1], "generating_index_fastq_from_s3")

        # Split by number of reads per worker (each read is composed of 4 lines)
        assert (num_lines % 4) == 0, "fastq file total number of lines is not multiple of 4!"
        num_reads = num_lines // 4
        reads_batch = ceil(num_reads / pipeline_params.fastq_chunks)
        read_pairs = [(reads_batch * i, (reads_batch * i) + reads_batch) for i in range(pipeline_params.fastq_chunks)]

        # Convert read pairs back to line numbers (starting in 1)
        line_pairs = [((l0 * 4) + 1, (l1 * 4) + 1) for l0, l1 in read_pairs]

        # Adjust last pair for num batches not multiple of number of total reads (last batch will have fewer lines)
        if line_pairs[-1][1] > num_lines:
            l0, _ = line_pairs[-1]
            line_pairs[-1] = (l0, num_lines + 1)

        # Get byte ranges from line pairs using GZip index
        # TODO maybe call this function using lithops? currenty it needs to download the tab file to local process, compare if it is faster invoking remote lambda
        byte_ranges = get_ranges_from_line_pairs(pipeline_params, lithops, line_pairs, subStat)
        chunks = [
            {"chunk_id": i, "line_0": line_0, "line_1": line_1, "range_0": range_0, "range_1": range_1}
            for i, ((line_0, line_1), (range_0, range_1)) in enumerate(zip(line_pairs, byte_ranges))
        ]
        subStat.timer_stop("prepare_fastq_chunks")
        return chunks, subStat
    elif pipeline_params.fastq_sra is not None:
        # fastq-dump works by number of reads, number of lines = number of reads * 4
        num_reads = get_sra_metadata(pipeline_params)
        reads_batch = ceil(num_reads / pipeline_params.fastq_chunks)
        read_pairs = [
            (reads_batch * i + (1 if i == 0 else 0), (reads_batch * i) + reads_batch)
            for i in range(pipeline_params.fastq_chunks)
        ]

        # Adjust last pair for num batches not multiple of number of total reads (last batch will have fewer reads)
        if read_pairs[-1][1] > num_reads:
            l0, _ = read_pairs[-1]
            read_pairs[-1] = (l0, num_reads)

        chunks = [{"chunk_id": i, "read_0": read_0, "read_1": read_1} for i, (read_0, read_1) in enumerate(read_pairs)]
        subStat.timer_stop("prepare_fastq_chunks")

        return chunks, subStat
    else:
        raise Exception("fastq reference required")


def get_sra_metadata(pipeline_params: PipelineRun) -> int:
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {"db": "sra", "id": pipeline_params.fastq_sra, "retmode": "xml"}

    response = requests.get(base_url, params=params)

    if response.status_code == 200:
        xml_data = response.text
        root = ET.fromstring(xml_data)

        for run in root.iter("RUN"):
            reads = int(run.get("total_spots"))
            logger.debug("Read total reads from efetch for sequence %s ", pipeline_params.fastq_sra)
            return reads
    else:
        raise Exception(f"Error fetching metadata for {pipeline_params.fastq_sra}: {response.status_code}")
