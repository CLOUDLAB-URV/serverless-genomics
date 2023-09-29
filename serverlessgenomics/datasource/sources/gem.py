from __future__ import annotations

import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ...pipeline import PipelineParameters


def get_gem_chunk_storage_key(pipeline_params: PipelineParameters, fasta_chunk_id: int) -> str:
    return os.path.join(
        pipeline_params.gem_index_prefix,
        pipeline_params.fasta_path.key,
        f"{pipeline_params.fasta_chunks}-chunks",
        "chunk" + str(fasta_chunk_id).zfill(4) + ".gem",
    )


def get_gem_chunk_storage_prefix(pipeline_params: PipelineParameters) -> str:
    return os.path.join(
        pipeline_params.gem_index_prefix,
        pipeline_params.fasta_path.key,
        f"{pipeline_params.fasta_chunks}-chunks",
    )
