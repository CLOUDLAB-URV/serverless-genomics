from __future__ import annotations

import io
import logging
import re
import subprocess
import os
import tempfile
import time
from math import ceil
from typing import TYPE_CHECKING, List, Tuple

import numpy as np
import pandas as pd

from ..utils import force_delete_path, S3Path, try_head_object

if TYPE_CHECKING:
    from ..parameters import PipelineRun, Lithops
    from mypy_boto3_s3 import S3Client
    import lithops

logger = logging.getLogger(__name__)

CHUNK_SIZE = 65536
RE_WINDOWS = re.compile(r'#\d+: @ \d+ / \d+ L\d+ \( \d+ @\d+ \)')
RE_NUMS = re.compile(r'\d+')
RE_NLINES = re.compile(r'Number of lines\s+:\s+\d+')


def _get_gztool_path():
    """
    Utility function that returns the absolute path for gzip file binary or raises exception if it is not found
    """
    proc = subprocess.run(['which', 'gztool'], check=True, capture_output=True, text=True)
    path = proc.stdout.rstrip('\n')
    logger.debug('Using gztool located in %s', path)
    return path


# @lithops_callee
def generate_idx_from_gzip(pipeline_params: PipelineRun, gzip_file_path: S3Path, storage: lithops.Storage):
    """
    Lithops callee function
    Create index file from gzip archive using gztool (https://github.com/circulosmeos/gztool)
    """
    gztool = _get_gztool_path()
    tmp_index_file_name = tempfile.mktemp()
    s3: S3Client = storage.get_client()

    try:
        res = s3.get_object(Bucket=gzip_file_path.bucket, Key=gzip_file_path.key)
        data_stream = res['Body']
        force_delete_path(tmp_index_file_name)
        t0 = time.perf_counter()

        # Create index and save to tmp file
        index_proc = subprocess.Popen([gztool, '-i', '-x', '-I', tmp_index_file_name],
                                      stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        chunk = data_stream.read(CHUNK_SIZE)
        while chunk != b"":
            index_proc.stdin.write(chunk)
            chunk = data_stream.read(CHUNK_SIZE)
        if hasattr(data_stream, 'close'):
            data_stream.close()

        stdout, stderr = index_proc.communicate()
        # logger.debug(stdout.decode('utf-8'))
        # logger.debug(stderr.decode('utf-8'))
        if index_proc.returncode > 0:
            logger.debug(stdout.decode('utf-8'))
            logger.debug(stderr.decode('utf-8'))
            raise Exception('Error creating gz index')

        # Generate list of access windows from index
        proc = subprocess.run([gztool, '-ell', '-I', tmp_index_file_name], check=True,
                              capture_output=True, text=True)
        output = proc.stdout

        # Store index binary file
        index_key, tab_key = pipeline_params.fastqgz_idx_keys
        s3.upload_file(Filename=tmp_index_file_name, Bucket=pipeline_params.storage_bucket, Key=index_key)

        # Get the total number of lines
        total_lines: int = int(RE_NUMS.findall(RE_NLINES.findall(output).pop()).pop())
        logger.debug('Indexed gzipped text file with %d total lines', total_lines)
        t1 = time.perf_counter()
        logger.debug('Index generated in %.3f seconds', t1 - t0)

        # Generator function that parses output to avoid copying all window data as lists
        def _lines_generator():
            for f in RE_WINDOWS.finditer(output):
                # nums = np.array([int(n) for n in RE_NUMS.findall(f.group())], dtype=np.int)
                nums = [int(n) for n in RE_NUMS.findall(f.group())]
                yield nums

        # Generate data frame that stores gzip index windows offsets

        df = pd.DataFrame(_lines_generator(),
                          columns=['window', 'compressed_byte', 'uncompressed_byte',
                                   'line_number', 'window_size', 'window_offset'])
        df.set_index(['window'], inplace=True)

        # Store data frame as parquet
        out_stream = io.BytesIO()
        df.to_parquet(out_stream, engine='pyarrow')
        out_stream.seek(0)

        s3.put_object(Bucket=pipeline_params.storage_bucket, Key=tab_key, Body=out_stream,
                      Metadata={'total_lines': str(total_lines)})
        return total_lines
    finally:
        force_delete_path(tmp_index_file_name)


def get_ranges_from_line_pairs(pipeline_params: PipelineRun, lithops: Lithops, pairs: List[Tuple[int, int]]):
    _, tab_key = pipeline_params.fastqgz_idx_keys
    meta_obj_body = lithops.storage.get_object(bucket=pipeline_params.storage_bucket, key=tab_key)
    meta_buff = io.BytesIO(meta_obj_body)
    meta_buff.seek(0)
    df = pd.read_parquet(meta_buff)
    line_indexes = df['line_number'].to_numpy()
    num_windows = df.shape[0]

    head_fastq = lithops.storage.head_object(bucket=pipeline_params.fastq_path.bucket,
                                             key=pipeline_params.fastq_path.key)
    fastqgz_sz = head_fastq['content-length']

    byte_ranges = [None] * len(pairs)
    for i, (line_0, line_1) in enumerate(pairs):
        # Find the closest window index for line_0
        window_head_idx = (np.abs(line_indexes - line_0)).argmin()
        # Check if window line entry pont is past requested line_0, if so, get previous window
        window_head_line = df.iloc[window_head_idx]['line_number']
        if window_head_line > line_0:
            window_head_idx = window_head_idx - 1
        # Get offset in compressed archive for window 0
        window0_offset = df.iloc[window_head_idx]['compressed_byte']

        # Find the closest window index for line_0
        widow_tail_idx = (np.abs(line_indexes - line_1)).argmin()
        # Check if window line entry pont is before requested line_1, if so, get next window
        window_tail_line = df.iloc[widow_tail_idx]['line_number']
        if window_tail_line < line_1:
            widow_tail_idx = widow_tail_idx + 1
        if widow_tail_idx >= num_windows:
            # Adjust offset for lines inside last window, use end of compressed archive for 2nd offset
            window1_offset = fastqgz_sz
        else:
            window1_offset = df.iloc[widow_tail_idx]['compressed_byte']

        byte_ranges[i] = (window0_offset, window1_offset)

    return byte_ranges


def generate_fastqgz_index_from_s3(pipeline_params: PipelineRun, lithops: Lithops) -> int:
    """
    Generate gzip index file for FASTQ in storage if necessary, returns number of lines of FASTQ file
    """
    # Check if fastq file exists
    fastq_head = try_head_object(lithops.storage, pipeline_params.fastq_path.bucket, pipeline_params.fastq_path.key)
    if fastq_head is None:
        raise Exception(f'fastq file with key {pipeline_params.fastq_path} does not exists')

    # Check if fastqgz index file exists
    index_key, tab_key = pipeline_params.fastqgz_idx_keys
    fastq_idx_head = try_head_object(lithops.storage,
                                     pipeline_params.fastq_path.bucket, index_key)
    fastq_tab_head = try_head_object(lithops.storage,
                                     pipeline_params.fastq_path.bucket, tab_key)
    if fastq_idx_head is None or fastq_tab_head is None:
        # Generate gzip index file for compressed fastq input
        logger.info('Generating gzip index file for FASTQ %s', pipeline_params.fastq_path.stem)
        total_lines = lithops.invoker.call('generate_fasta_idx', pipeline_params,
                                           generate_idx_from_gzip, (pipeline_params, pipeline_params.fastq_path))
    else:
        # Get total lines from header metadata
        logger.debug('Fastqz index for %s found', pipeline_params.fastq_path.stem)
        total_lines = int(fastq_tab_head['x-amz-meta-total_lines'])
    logger.info('Read %d sequences from FASTQ %s', total_lines / 4, pipeline_params.fastq_path.stem)

    return total_lines


def prepare_fastq_chunks(pipeline_params: PipelineRun, lithops: Lithops):
    """
    Calculate fastq byte ranges for chunks of a pipeline run, generate gzip index if needed
    """
    if pipeline_params.fastq_path is not None:
        num_lines = generate_fastqgz_index_from_s3(pipeline_params, lithops)
    elif pipeline_params.fastq_sra is not None:
        # TODO implement get fastq file from sra archive
        raise NotImplementedError()
    else:
        raise Exception('fastq reference required')

    # Split by number of reads per worker (each read is composed of 4 lines)
    assert (num_lines % 4) == 0, 'fastq file total number of lines is not multiple of 4!'
    num_reads = num_lines // 4
    reads_batch = ceil(num_reads / pipeline_params.fastq_chunks)
    read_pairs = [(reads_batch * i, (reads_batch * i) + reads_batch) for i in range(pipeline_params.fastq_chunks)]

    # Convert read pairs back to line numbers (starting in 1)
    line_pairs = [((l0 * 4) + 1, (l1 * 4) + 1) for l0, l1 in read_pairs]

    # Adjust last pair for num batches not multiple of number of total reads (last batch will have fewer lines)
    if line_pairs[-1][1] > num_lines:
        l0, _ = line_pairs[-1]
        line_pairs[-1] = (l0, num_lines)

    # Get byte ranges from line pairs using GZip index
    # TODO maybe call this function using lithops? currenty it needs to download the tab file to local process, compare if it is faster invoking remote lambda
    byte_ranges = get_ranges_from_line_pairs(pipeline_params, lithops, line_pairs)
    return byte_ranges


# TODO remake function, not optimal...
def fastq_to_mapfun(fastq_file_key: str, fastq_chunk_data: str) -> str:
    """
    Function executed within the map function to retrieve the relevant fastq chunk from object storage
    """
    seq_name = fastq_file_key

    subprocess.call(['chmod', '+x', 'fastq-dump'])
    subprocess.run(['vdb-config', '-i'])  # To supress a warning that appears the first time vdb-config is used

    # Report cloud identity so it can take data from s3 needed to be executed only once per vm
    output = str(subprocess.run(['vdb-config', '--report-cloud-identity', 'yes'], capture_output=True).stdout)

    os.chdir(f"/tmp")
    temp_fastq = f'/tmp/' + seq_name + f'_chunk{fastq_chunk_data["number"]}.fastq'
    data_output = subprocess.run(['fastq-dump', str(seq_name), '-X', str(int(fastq_chunk_data["start_line"])), '-N',
                                  str(int(fastq_chunk_data["end_line"])), '-O', f'/tmp'],
                                 capture_output=True)
    os.rename(f'/tmp/' + seq_name + '.fastq', temp_fastq)

    return temp_fastq
