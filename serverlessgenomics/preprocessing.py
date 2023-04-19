from __future__ import annotations

import logging
from math import ceil
from typing import TYPE_CHECKING

from serverlessgenomics.datasource.sources.fastqgz import check_fastqgz_index, get_ranges_from_line_pairs
from serverlessgenomics.datasource.sources.sra import get_sra_metadata
from .datasource.sources.fasta import generate_faidx_from_s3, get_fasta_byte_ranges
from .stats import Stats

if TYPE_CHECKING:
    from .pipelineparams import PipelineParameters, Lithops

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
        chunks = [
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
        logger.info("Generated %d chunks for %s", len(chunks), pipeline_params.fastq_path)
        subStat.timer_stop("prepare_fastq_chunks")
        return chunks, subStat
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

        chunks = [
            {"source": "sra", "chunk_id": i, "read_0": read_0, "read_1": read_1}
            for i, (read_0, read_1) in enumerate(read_pairs)
        ]
        subStat.timer_stop("prepare_fastq_chunks")

        return chunks, subStat
    else:
        raise Exception("fastq reference required")


def prepare_fasta_chunks(pipeline_params: PipelineParameters, lithops: Lithops):
    """
    Calculate fasta byte ranges and metadata for chunks of a pipeline run, generate faidx index if needed
    """
    subStat = Stats()
    # Get number of sequences from fasta file, generate faidx file if needed
    subStat.timer_start("prepare_fasta_chunks")
    num_sequences = generate_faidx_from_s3(pipeline_params, lithops, subStat)
    fasta_chunks = get_fasta_byte_ranges(pipeline_params, lithops, num_sequences)
    subStat.timer_stop("prepare_fasta_chunks")

    return fasta_chunks, subStat