from __future__ import annotations

import logging
import multiprocessing
import os
import subprocess
import tempfile
from math import ceil
from time import time
from typing import TYPE_CHECKING

from lithops.storage.utils import StorageNoSuchKeyError

from serverlessgenomics.datasource.sources.fastqgz import check_fastqgz_index, get_ranges_from_line_pairs
from serverlessgenomics.datasource.sources.sra import get_sra_metadata
from .datasource import fetch_fasta_chunk
from .datasource.sources.fasta import generate_faidx_from_s3, get_fasta_byte_ranges
from .stats import Stats
from .utils import force_delete_local_path

if TYPE_CHECKING:
    from typing import List
    from .pipeline import PipelineParameters, Lithops
    from lithops import Storage


logger = logging.getLogger(__name__)


def prepare_fastq_chunks(pipeline_params: PipelineParameters, lithops: Lithops):
    """
    Generate FASTQ chunks metadata
    """
    subStat = Stats()
    subStat.timer_start("prepare_fastq_chunks")
    if pipeline_params.fastq_path is not None:
        # FASTQ source is S3
        num_lines = check_fastqgz_index(pipeline_params, lithops)
        logger.info("Read %d sequences from %s", num_lines / 4, pipeline_params.fastq_path)

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
        logger.info(
            "Calculating byte ranges of %s for %d chunks...", pipeline_params.fastq_path, pipeline_params.fastq_chunks
        )
        byte_ranges = lithops.invoker.call(get_ranges_from_line_pairs, (pipeline_params, line_pairs))
        fastq_chunks = [
            {
                "source": "s3_fastqgzip",
                "chunk_id": i,
                "line_0": line_0,
                "line_1": line_1,
                "range_0": range_0,
                "range_1": range_1,
            }
            for i, ((line_0, line_1), (range_0, range_1)) in enumerate(zip(line_pairs, byte_ranges))
        ]
        logger.info("Generated %d chunks for %s", len(fastq_chunks), pipeline_params.fastq_path)
        subStat.timer_stop("prepare_fastq_chunks")
    elif pipeline_params.sra_accession is not None:
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

        fastq_chunks = [
            {"source": "sra", "chunk_id": i, "read_0": read_0, "read_1": read_1}
            for i, (read_0, read_1) in enumerate(read_pairs)
        ]
        subStat.timer_stop("prepare_fastq_chunks")
    else:
        raise Exception("fastq reference required")

    if pipeline_params.fastq_chunk_range is not None:
        # Compute only specified FASTQ chunk range
        r0, r1 = pipeline_params.fastq_chunk_range
        logger.info("Using only FASTQ chunks in range %s", pipeline_params.fastq_chunk_range.__repr__())
        fastq_chunks = fastq_chunks[r0:r1]

    return fastq_chunks, subStat


def prepare_fasta_chunks(pipeline_params: PipelineParameters, lithops: Lithops):
    """
    Calculate fasta byte ranges and metadata for chunks of a pipeline run, generate faidx index if needed
    """
    subStat = Stats()
    # Get number of sequences from fasta file, generate faidx file if needed
    subStat.timer_start("prepare_fasta_chunks")
    num_sequences = generate_faidx_from_s3(pipeline_params, lithops, subStat)
    fasta_chunks = get_fasta_byte_ranges(pipeline_params, lithops, num_sequences)
    if pipeline_params.fasta_chunk_range is not None:
        # Compute only specified FASTA chunk range
        r0, r1 = pipeline_params.fasta_chunk_range
        logger.info("Using only FASTA chunks in range %s", pipeline_params.fasta_chunk_range.__repr__())
        fasta_chunks = fasta_chunks[r0:r1]
    logger.info("Generated %d chunks for %s", len(fasta_chunks), pipeline_params.fasta_path.as_uri())
    subStat.timer_stop("prepare_fasta_chunks")

    return fasta_chunks, subStat


def generate_gem_indexer_iterdata(
    pipeline_params: PipelineParameters, run_id: str, fasta_chunks: List[dict]
) -> List[dict]:
    iterdata = []

    for fa_i, fa_ch in enumerate(fasta_chunks):
        params = {"pipeline_params": pipeline_params, "run_id": run_id, "fasta_chunk_id": fa_i, "fasta_chunk": fa_ch}
        iterdata.append(params)

    return iterdata


def gem_indexer(
    pipeline_params: PipelineParameters, run_id: str, fasta_chunk_id: int, fasta_chunk: dict, storage: Storage
):
    # Stats
    stat, timestamps, data_size = Stats(), Stats(), Stats()
    stat.timer_start(fasta_chunk_id)
    timestamps.store_size_data("start", time())

    # Initialize names
    gem_index_filename = os.path.join(f"{fasta_chunk_id}.gem")

    # Check if gem file already exists
    try:
        storage.head_object(bucket=pipeline_params.storage_bucket, key=f"gem/{gem_index_filename}")
        stat.timer_stop(fasta_chunk_id)
        return stat.get_stats()
    except StorageNoSuchKeyError:
        pass

    # Make temp dir and ch into it, save pwd to restore it later
    tmp_dir = tempfile.mkdtemp()
    pwd = os.getcwd()
    os.chdir(tmp_dir)

    try:
        # Get fasta chunk and store it to disk in tmp directory
        timestamps.store_size_data("download_fasta", time())
        fasta_chunk_filename = f"chunk_{fasta_chunk['chunk_id']}.fasta"
        fetch_fasta_chunk(fasta_chunk, fasta_chunk_filename, storage, pipeline_params.fasta_path)
        data_size.store_size_data(fasta_chunk_filename, os.path.getsize(fasta_chunk_filename) / (1024 * 1024))

        # gem-indexer appends .gem to output file
        timestamps.store_size_data("gem_indexer", time())
        cmd = [
            "gem-indexer",
            "--input",
            fasta_chunk_filename,
            "--threads",
            str(multiprocessing.cpu_count()),
            "-o",
            gem_index_filename.replace(".gem", ""),
        ]
        print(" ".join(cmd))
        out = subprocess.run(cmd, capture_output=True)
        print(out.stderr.decode("utf-8"))

        # TODO apparently 1 is return code for success (why)
        assert out.returncode == 1

        # Upload the gem file to storage
        timestamps.store_size_data("upload_gem", time())
        storage.upload_file(
            file_name=gem_index_filename, bucket=pipeline_params.storage_bucket, key=f"gem/{gem_index_filename}"
        )

        # Finish stats
        data_size.store_size_data(gem_index_filename, os.path.getsize(gem_index_filename) / (1024 * 1024))
        stat.timer_stop(fasta_chunk_id)
        timestamps.store_size_data("end", time())

        # Store dict
        stat.store_dictio(timestamps.get_stats(), "timestamps", fasta_chunk_id)
        stat.store_dictio(data_size.get_stats(), "data_sizes", fasta_chunk_id)
        return stat.get_stats()
    finally:
        os.chdir(pwd)
        force_delete_local_path(tmp_dir)
