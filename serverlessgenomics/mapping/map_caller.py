from __future__ import annotations

import collections
import logging
from typing import TYPE_CHECKING

from .alignment_mapper import align_mapper, index_correction, filtered_index_to_mpileup
from ..pipeline import PipelineParameters, Lithops, PipelineRun

if TYPE_CHECKING:
    from typing import Tuple

logger = logging.getLogger(__name__)


def format_align_mapper_id(fasta_chunk_id: int, fastq_chunk_id: int) -> str:
    return "fa" + str(fasta_chunk_id).zfill(4) + "-" + "fq" + str(fastq_chunk_id).zfill(4)


def unformat_align_mapper_id(mapper_id: str) -> Tuple[int, int]:
    """
    Returns tuple (fasta_chunk_id, fastq_chunk_id) from mapper formatted id
    """
    fa, fq = mapper_id.split("-")
    return int(fa.replace("fa", "")), int(fq.replace("fq", ""))


def format_index_correction_mapper_id(fastq_chunk_id: int) -> str:
    return "fq" + str(fastq_chunk_id).zfill(4)


def unformat_index_correction_mapper_id(index_correction_mapper_id: str) -> int:
    return int(index_correction_mapper_id.replace("fq", ""))


def generate_align_mapping_iterdata(pipeline_params: PipelineParameters, pipeline_run: PipelineRun):
    iterdata = [
        {
            "pipeline_params": pipeline_params,
            "run_id": pipeline_run.run_id,
            "mapper_id": format_align_mapper_id(fa_ch["chunk_id"], fq_ch["chunk_id"]),
            "fasta_chunk": fa_ch,
            "fastq_chunk": fq_ch,
        }
        for fq_ch in pipeline_run.fastq_chunks
        for fa_ch in pipeline_run.fasta_chunks
    ]

    return iterdata


def generate_index_correction_iterdata(pipeline_params, pipeline_run):
    # Group gem mapper output by fastq chunk id
    grouped_fastq_mappers = collections.defaultdict(list)
    for mapper_id, (map_key, _) in pipeline_run.alignment_maps.items():
        _, fastq_chunk_id = unformat_align_mapper_id(mapper_id)
        format_index_correction_mapper_id(fastq_chunk_id)
        grouped_fastq_mappers[fastq_chunk_id].append(map_key)

    iterdata = []

    for fq_id, map_keys in grouped_fastq_mappers.items():
        params = {
            "pipeline_params": pipeline_params,
            "run_id": pipeline_run.run_id,
            "mapper_id": format_index_correction_mapper_id(fq_id),
            "map_index_keys": map_keys,
        }
        iterdata.append(params)

    return iterdata


def generate_index_to_mpileup_iterdata(pipeline_params, pipeline_run):
    iterdata = []

    for fq_ch in pipeline_run.fastq_chunks:
        corrected_index_key = pipeline_run.corrected_indexes[format_index_correction_mapper_id(fq_ch["chunk_id"])]
        for fa_ch in pipeline_run.fasta_chunks:
            mapper_id = format_align_mapper_id(fa_ch["chunk_id"], fq_ch["chunk_id"])
            _, filtered_map_key = pipeline_run.alignment_maps[mapper_id]
            params = {
                "pipeline_params": pipeline_params,
                "run_id": pipeline_run.run_id,
                "mapper_id": mapper_id,
                "fasta_chunk": fa_ch,
                "filtered_map_key": filtered_map_key,
                "corrected_index_key": corrected_index_key,
            }
            iterdata.append(params)

    return iterdata


def run_full_alignment(pipeline_params: PipelineParameters, pipeline_run: PipelineRun, lithops: Lithops):
    """
    Execute the map phase
    """
    # subStat = Stats()

    # MAP: Stage 1
    logger.debug("PROCESSING MAP: STAGE 1")
    # subStat.timer_start("aligner_indexer")
    iterdata = generate_align_mapping_iterdata(pipeline_params, pipeline_run)
    align_mapper_result = lithops.invoker.map(align_mapper, iterdata)
    pipeline_run.alignment_maps = {
        mapper_id: (map_index_key, filtered_map_key)
        for mapper_id, map_index_key, filtered_map_key in align_mapper_result
    }

    # subStat.timer_stop("aligner_indexer")
    # align_mapper_result, timers = split_data_result(align_mapper_result)
    # subStat.store_dictio(timers, "function_details", "aligner_indexer")

    # MAP: Index correction
    logger.debug("PROCESSING INDEX CORRECTION")
    # subStat.timer_start("index_correction")

    iterdata = generate_index_correction_iterdata(pipeline_params, pipeline_run)
    index_correction_result = lithops.invoker.map(index_correction, iterdata)
    pipeline_run.corrected_indexes = {
        mapper_id: corrected_index_key for mapper_id, corrected_index_key in index_correction_result
    }

    # subStat.timer_stop("index_correction")
    # index_correction_result, timers = split_data_result(index_correction_result)
    # subStat.store_dictio(timers, "function_details", "index_correction")

    # Map: Stage 2
    logger.debug("PROCESSING MAP: STAGE 2")
    # subStat.timer_start("filter_index_to_mpileup")
    iterdata = generate_index_to_mpileup_iterdata(pipeline_params, pipeline_run)
    index_to_mpileup_result = lithops.invoker.map(filtered_index_to_mpileup, iterdata)
    pipeline_run.aligned_mpileups = {mapper_id: mpileup_key for mapper_id, mpileup_key in index_to_mpileup_result}

    # subStat.timer_stop("filter_index_to_mpileup")
    # alignment_output, timers = split_data_result(alignment_output)
    # subStat.store_dictio(timers, "function_details", "filter_index_to_mpileup")
