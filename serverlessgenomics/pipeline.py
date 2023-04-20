from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional, Tuple

import lithops
import uuid

from .lithopswrapper import LithopsInvokerWrapper
from .utils import S3Path, guess_sra_accession_from_fastq_path, validate_sra_accession_id

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PipelineParameters:
    """
    Dataclass to store a pipeline's input parameters and configuration
    """

    # FASTA parameters (reference genome)
    # Storage path for FASTA (reference genome) input file
    fasta_path: S3Path
    # Number of chunks to split FASTA input file into
    fasta_chunks: Optional[int] = None

    # FASTQ parameters (sequence read)
    # Storage path for fastq input file
    fastq_path: Optional[S3Path] = None
    # SRA run accession for FASTQ sequence read input
    sra_accession: Optional[str] = None
    # Number of chunks to split fastq input file into
    fastq_chunks: Optional[int] = None

    # Variant Calling parameters
    # TODO what is tolerance? (ask Lucio)
    tolerance: int = 0

    # Debug parameters
    # fastq chunks to be processed
    fastq_chunk_range: Tuple[int, int] = None
    # fasta chunks to be processed
    fasta_chunk_range: Tuple[int, int] = None
    # Skip PreProcessing Stage
    skip_prep: bool = False
    # Skip Map Stage
    skip_map: bool = False
    # Skip Reduce Stage
    skip_reduce: bool = False

    # Lithops settings
    lithops_settings: dict = None

    # Bucket name with write permissions to store preprocessed, intermediate and output data
    storage_bucket: str = "serverless-genomics"
    # Prefix for storing cached fastqgz indexes
    fastqgz_idx_prefix: str = "fastqgz-indexes/"
    # Prefix for storing cached faidx indexes
    faidx_prefix: str = "faidx-indexes/"
    # Prefix for storing cached reference genome indexes
    gem_index_prefix: str = "gem-indexes/"
    # Prefix for output data results
    output_prefix: str = "output/"

    # Log level
    log_level: str = "INFO"
    # Log stats
    log_stats: bool = False


@dataclass
class PipelineRun:
    """
    Dataclass to store a Pipeline execution state
    """

    # Run input parameters
    parameters: PipelineParameters
    # Run ID
    run_id: str

    fastq_chunks = None
    fasta_chunks = None
    gem_chunk_ids = None
    alignment_batches = None


@dataclass(frozen=True)
class Lithops:
    """
    Dataclass to encapsulate Lithops function executor and storage clients from a single session
    """

    storage: lithops.Storage
    invoker: LithopsInvokerWrapper


def validate_parameters(params: dict) -> PipelineParameters:
    """
    Validate and populate missing input parameters. Returns a correct PipelineParameters
    dataclass instance.
    """

    if "fasta_path" not in params:
        raise KeyError("fasta_path")
    if "fasta_chunks" not in params:
        raise KeyError("fasta_chunks")

    params["fasta_path"] = S3Path.from_uri(params["fasta_path"])

    if "fastq_path" in params:
        # Get fastq sequence read from S3
        params["fastq_path"] = S3Path.from_uri(params["fastq_path"])
        sra_accession = guess_sra_accession_from_fastq_path(params["fastq_path"].as_uri())

        if "sra_accession" in params and sra_accession is not None:
            # Check if guessed SRA accession ID matches
            assert params["sra_accession"] == sra_accession, "Specified SRA accession ID does not " "match ID in path"
        if "sra_accession" not in params:
            assert sra_accession is not None, "Could not guess SRA accession form fastq path, specify it explicitly"
            logging.info("Guessed %s as SRA accession id from fastq path", sra_accession)
            params["sra_accession"] = sra_accession
    else:
        # Get fastq read from SRA
        try:
            assert "sra_accession" in params
            assert validate_sra_accession_id(params["sra_accession"])
        except AssertionError as e:
            logging.error("FASTQ S3 path or valid SRA accession ID are required")
            raise e

    return PipelineParameters(**params)


def new_pipeline_run(pipeline_parameters: PipelineParameters) -> PipelineRun:
    run_id = str(uuid.uuid4())
    run = PipelineRun(parameters=pipeline_parameters, run_id=run_id)

    logger.info(
        "Created new run\n######################################################\n"
        "Pipeline Run ID = %s\n"
        "######################################################",
        run.run_id,
    )
    return run
