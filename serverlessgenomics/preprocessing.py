from __future__ import annotations

import collections
import logging
import multiprocessing
import os
import re
import subprocess
import tempfile
from math import ceil
from time import time
from typing import TYPE_CHECKING, Set

from lithops.storage.utils import StorageNoSuchKeyError

from .datasource.sources.fastqgz import check_fastqgz_index, get_ranges_from_line_pairs
from .datasource.sources.sra import get_sra_metadata
from .datasource.datasources import FASTQSource, FASTASource
from .datasource import fetch_fasta_chunk
from .datasource.sources.fasta import generate_faidx_from_s3, get_fasta_byte_ranges
from .datasource.sources.gem import get_gem_chunk_storage_key, get_gem_chunk_storage_prefix, gem_indexer
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
                "source": FASTQSource.S3_GZIP,
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
            {"source": FASTQSource.SRA, "chunk_id": i, "read_0": read_0, "read_1": read_1}
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


def prepare_gem_chunks(pipeline_params: PipelineParameters, fasta_chunks: list[dict], lithops: Lithops) -> Set[int]:
    """
    Generate GEM indexed file metadata
    """
    # Check if all GEM files exist for the input FASTA file and specified number of chunks
    gems_prefix = get_gem_chunk_storage_prefix(pipeline_params)

    cached_gems_keys = lithops.storage.list_keys(bucket=pipeline_params.storage_bucket, prefix=gems_prefix)
    cached_gem_chunk_ids = []

    for cached_gems_key in cached_gems_keys:
        basename = os.path.basename(cached_gems_key)
        print(basename)
        matches = re.findall(r'\d+', basename)
        assert len(matches) == 1
        chunk_id = int(matches.pop())
        cached_gem_chunk_ids.append(chunk_id)

    cached_gem_chunk_ids = set(cached_gem_chunk_ids)
    requested_gems_ids = set(range(pipeline_params.fasta_chunks))

    if cached_gem_chunk_ids == requested_gems_ids:
        # All chunks are already in storage
        logger.info("Using %d cached GEM files in storage (prefix=\"%s\")", len(cached_gem_chunk_ids), gems_prefix)
        return cached_gem_chunk_ids

    if not cached_gem_chunk_ids:
        # All chunks are missing
        iterdata = generate_gem_indexer_iterdata(pipeline_params, fasta_chunks)
        logger.debug("All GEM chunks are missing (%d in total)", len(requested_gems_ids))
    else:
        # Only some chunks are missing
        missing_chunks_ids = requested_gems_ids - cached_gem_chunk_ids
        logger.debug("Some GEM chunks are missing (%s)", missing_chunks_ids.__repr__())
        missing_fasta_chunks = []
        for fa_ch in fasta_chunks:
            if fa_ch["chunk_id"] not in missing_chunks_ids:
                missing_fasta_chunks.append(fa_ch)
        iterdata = generate_gem_indexer_iterdata(pipeline_params, missing_fasta_chunks)

    logger.info("Going to index %d GEM chunks", len(iterdata))
    lithops.invoker.map(gem_indexer, iterdata)

    return requested_gems_ids


def generate_gem_indexer_iterdata(pipeline_params: PipelineParameters, fasta_chunks: List[dict]) -> List[dict]:
    iterdata = []

    for fa_ch in fasta_chunks:
        params = {"pipeline_params": pipeline_params, "fasta_chunk_id": fa_ch["chunk_id"], "fasta_chunk": fa_ch}
        iterdata.append(params)

    return iterdata
