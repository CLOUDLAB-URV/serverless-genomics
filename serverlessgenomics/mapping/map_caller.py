import collections
import copy
import logging

from .alignment_mapper import align_mapper, index_correction, filter_index_to_mpileup, gem_indexer
from ..pipeline import PipelineParameters, Lithops, PipelineRun
from ..stats import Stats
from ..utils import split_data_result

logger = logging.getLogger(__name__)


def generate_gem_indexer_iterdata(pipeline_params: PipelineParameters, fasta_chunks):
    iterdata = []

    for fa_i, fa_ch in enumerate(fasta_chunks):
        params = {"pipeline_params": pipeline_params, "fasta_chunk_id": fa_i, "fasta_chunk": fa_ch}
        iterdata.append(params)

    return iterdata


def generate_align_mapping_iterdata(pipeline_params: PipelineParameters, fasta_chunks, fastq_chunks):
    iterdata = []
    for fq_i, fq_ch in enumerate(fastq_chunks):
        for fa_i, fa_ch in enumerate(fasta_chunks):
            params = {
                "pipeline_params": pipeline_params,
                "fasta_chunk_id": fa_i,
                "fasta_chunk": fa_ch,
                "fastq_chunk": fq_ch,
                "fastq_chunk_id": fq_i,
            }
            iterdata.append(params)
    return iterdata


def generate_index_correction_iterdata(pipeline_params, gem_mapper_output):
    # Group gem mapper output by fastq chunk id
    fq_groups = collections.defaultdict(list)
    for fq, _, map_key, _ in gem_mapper_output:
        fq_groups[fq].append(map_key)

    iterdata = [
        {"pipeline_params": pipeline_params, "fastq_chunk_id": fq_id, "map_index_keys": map_keys}
        for fq_id, map_keys in fq_groups.items()
    ]

    return iterdata


def generate_index_to_mpileup_iterdata(
    pipeline_params, fasta_chunks, fastq_chunks, gem_mapper_output, corrected_indexes
):
    iterdata = []

    # Convert corrected index output (list of tuples) to sorted list by fastq chunk id
    corrected_indexes_fq = [tup[1] for tup in sorted(corrected_indexes, key=lambda tup: tup[0])]

    for fq_i, fa_i, _, filter_map_index in gem_mapper_output:
        iterdata.append(
            {
                "pipeline_params": pipeline_params,
                "fasta_chunk_id": fa_i,
                "fasta_chunk": fasta_chunks[fa_i],
                "fastq_chunk_id": fq_i,
                "fastq_chunk": fastq_chunks[fq_i],
                "filtered_map_key": filter_map_index,
                "corrected_index_key": corrected_indexes_fq[fq_i],
            }
        )

    return iterdata


def run_full_alignment(pipeline_params: PipelineParameters, pipeline_run: PipelineRun, lithops: Lithops):
    """
    Execute the map phase
    """
    subStat = Stats()

    # MAP: Stage 0.1
    logger.debug("PROCESSING GEM")
    iterdata = generate_gem_indexer_iterdata(pipeline_params, pipeline_run.fasta_chunks)
    subStat.timer_start("gem_generator")
    timers = lithops.invoker.map(gem_indexer, iterdata)
    subStat.timer_stop("gem_generator")
    subStat.store_dictio(timers, "function_details", "gem_generator")

    # MAP: Stage 1
    logger.debug("PROCESSING MAP: STAGE 1")
    subStat.timer_start("aligner_indexer")
    iterdata = generate_align_mapping_iterdata(pipeline_params, pipeline_run.fasta_chunks, pipeline_run.fastq_chunks)
    aligner_indexer_result = lithops.invoker.map(align_mapper, iterdata)
    subStat.timer_stop("aligner_indexer")
    aligner_indexer_result, timers = split_data_result(aligner_indexer_result)
    subStat.store_dictio(timers, "function_details", "aligner_indexer")

    # MAP: Index correction
    logger.debug("PROCESSING INDEX CORRECTION")
    subStat.timer_start("index_correction")

    iterdata = generate_index_correction_iterdata(pipeline_params, aligner_indexer_result)
    index_correction_result = lithops.invoker.map(index_correction, iterdata)
    subStat.timer_stop("index_correction")
    index_correction_result, timers = split_data_result(index_correction_result)
    subStat.store_dictio(timers, "function_details", "index_correction")

    # Map: Stage 2
    logger.debug("PROCESSING MAP: STAGE 2")
    subStat.timer_start("filter_index_to_mpileup")
    iterdata = generate_index_to_mpileup_iterdata(
        pipeline_params,
        pipeline_run.fasta_chunks,
        pipeline_run.fastq_chunks,
        aligner_indexer_result,
        index_correction_result,
    )
    alignment_output = lithops.invoker.map(filter_index_to_mpileup, iterdata)
    subStat.timer_stop("filter_index_to_mpileup")
    alignment_output, timers = split_data_result(alignment_output)
    subStat.store_dictio(timers, "function_details", "filter_index_to_mpileup")

    return alignment_output, subStat
