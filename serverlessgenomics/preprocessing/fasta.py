from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from ..datasource.sources.fasta import generate_faidx_from_s3, get_fasta_byte_ranges

if TYPE_CHECKING:
    from ..pipeline import PipelineParameters, Lithops

logger = logging.getLogger(__name__)


def prepare_fasta_chunks(pipeline_params: PipelineParameters, lithops: Lithops):
    """
    Calculate fasta byte ranges and metadata for chunks of a pipeline run, generate faidx index if needed
    """
    # Get number of sequences from fasta file, generate faidx file if needed
    num_sequences = generate_faidx_from_s3(pipeline_params, lithops)
    fasta_chunks = get_fasta_byte_ranges(pipeline_params, lithops, num_sequences)

    if pipeline_params.fasta_chunk_range is not None:
        # Compute only specified FASTA chunk range
        r0, r1 = pipeline_params.fasta_chunk_range
        logger.info(
            "Using only FASTA chunks in range %s",
            pipeline_params.fasta_chunk_range.__repr__(),
        )
        fasta_chunks = fasta_chunks[r0:r1]

    logger.info("Generated %d chunks for %s", len(fasta_chunks), pipeline_params.fasta_path.as_uri())

    return fasta_chunks
