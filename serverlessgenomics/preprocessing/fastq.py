from __future__ import annotations

import logging
from math import ceil
from typing import TYPE_CHECKING, Set

from ..datasource.sources.fastqgz import check_fastqgz_index, get_ranges_from_line_pairs
from ..datasource.sources.sra import get_sra_metadata
from ..datasource.datasources import FASTQSource

if TYPE_CHECKING:
    from ..pipeline import PipelineParameters, Lithops

logger = logging.getLogger(__name__)


def prepare_fastq_chunks(pipeline_params: PipelineParameters, lithops: Lithops):
    """
    Generate FASTQ chunks metadata
    """
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
            "Calculating byte ranges of %s for %d chunks...",
            pipeline_params.fastq_path,
            pipeline_params.fastq_chunks,
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
            {
                "source": FASTQSource.SRA,
                "chunk_id": i,
                "read_0": read_0,
                "read_1": read_1,
            }
            for i, (read_0, read_1) in enumerate(read_pairs)
        ]
    else:
        raise Exception("fastq reference required")

    if pipeline_params.fastq_chunk_range is not None:
        # Compute only specified FASTQ chunk range
        r0, r1 = pipeline_params.fastq_chunk_range
        logger.info(
            "Using only FASTQ chunks in range %s",
            pipeline_params.fastq_chunk_range.__repr__(),
        )
        fastq_chunks = fastq_chunks[r0:r1]

    return fastq_chunks
