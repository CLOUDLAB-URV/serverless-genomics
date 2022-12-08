import logging
import time
from .alignment_mapper import gem_indexer_mapper, filter_index_to_mpileup
from ..parameters import PipelineRun, Lithops
import os
from .. import aux_functions as af

logger = logging.getLogger(__name__)

map1_cachefile = 'lithops_map1_checkpoint'
correction_cachefile = 'lithops_correction_checkpoint'
map2_cachefile = 'lithops_map2_checkpoint'


def generate_gem_indexer_mapper_iterdata(pipeline_params, fasta_chunks, fastq_chunks):
    iterdata = []
    for fa_i, fa_ch in enumerate(fasta_chunks):
        for fq_i, fq_ch in enumerate(fastq_chunks):
            params = {'pipeline_params': pipeline_params, 'mapper_id': f'fa{fa_i}-fq{fq_i}',
                      'fasta_chunk': fa_ch, 'fastq_chunk': fq_ch}
            iterdata.append(params)
    return iterdata


def run_full_alignment(pipeline_params: PipelineRun, lithops: Lithops, fasta_chunks, fastq_chunks):
    """
    Execute the map phase

    Args:
        pipeline_params (PipelineRun): pipeline arguments
        alignment_batches (list): iterdata generated in the preprocessing stage
        map_func (AlignmentMapper): class containing the map functions
        num_chunks (int): number of corrections needed

    Returns:
        float: time taken to execute this phase
    """
    # Delete old files
    # logger.debug('Deleting previous mapper outputs...')
    # af.delete_files(lithops.storage, pipeline_params,
    #                 cloud_prefixes=[f'{pipeline_params.file_format}/{pipeline_params.execution_name}/'])

    # logger.debug("Running Map Phase... %d functions", len(alignment_batches))

    # Initizalize execution debug info
    # start = time.perf_counter()

    # MAP: STAGE 1
    logger.debug("PROCESSING MAP: STAGE 1")

    iterdata = generate_gem_indexer_mapper_iterdata(pipeline_params, fasta_chunks, fastq_chunks)
    gem_indexer_mapper_result = lithops.invoker.map(gem_indexer_mapper, iterdata)

    # Load futures from previous execution
    # map1_futures = af.load_cache(map1_cachefile, pipeline_params)

    # Execute first map if futures were not found
    # iter_data = generate_alignment_iterdata
    # if not map1_futures:
    #     map1_futures = fexec.map(map_func.gem_indexer_mapper, alignment_batches,
    #                              timeout=int(pipeline_params.func_timeout_map))

    # Get results either from the old futures or the new execution
    # first_map_results = fexec.get_result(fs=map1_futures)

    # Dump futures into file
    # af.dump_cache(map1_cachefile, map1_futures, pipeline_params)

    ###################################################################
    #### MAP: GENERATE CORRECTED INDEXES
    ###################################################################
    # logger.debug("PROCESSING INDEX CORRECTION")

    # Load futures from previous execution
    # correction_futures = af.load_cache(correction_cachefile, pipeline_params)

    # Execute correction if futures were not found
    # if (not correction_futures):
    #     # Generate the iterdata for index correction
    #     index_iterdata = []
    #     for i in range(num_chunks):
    #         index_iterdata.append({'setname': pipeline_params.fq_seqname + '_fq' + str(i + 1),
    #                                'bucket': str(pipeline_params.storage_bucket),
    #                                'exec_param': pipeline_params.execution_name})
    #
    #     # Execute index correction
    #     correction_futures = fexec.map(index_correction_map, index_iterdata,
    #                                    timeout=int(pipeline_params.func_timeout_map))
    #
    # # Get results either from the old futures or the new execution
    # fexec.get_result(fs=correction_futures)
    #
    # # Dump futures into file
    # af.dump_cache(correction_cachefile, correction_futures, pipeline_params)

    ###################################################################
    #### MAP: STAGE 2
    ###################################################################
    # print("PROCESSING MAP: STAGE 2")
    #
    # # Load futures from previous execution
    # map2_futures = af.load_cache(map2_cachefile, pipeline_params)

    # Execute correction if futures were not found
    # if not map2_futures:
    #     # Generate new iterdata
    #     newiterdata = []
    #     for worker in first_map_results:
    #         newiterdata.append({
    #             'fasta_chunk': worker[0],
    #             'fastq_chunk': worker[1],
    #             'corrected_map_index_file': worker[2].split("-")[0] + ".txt",
    #             'filtered_map_file': worker[3],
    #             'base_name': worker[4],
    #             'old_id': worker[5],
    #             'exec_param': pipeline_params.execution_name
    #         })
    #
    #     # Execute second stage of map
    #     map2_futures = fexec.map(map_func.filter_index_to_mpileup, newiterdata,
    #                              timeout=int(pipeline_params.func_timeout_map))
    #
    # # Get results either from the old futures or the new execution
    # fexec.get_result(fs=map2_futures)
    #
    # # Dump futures into file
    # af.dump_cache(map2_cachefile, map2_futures, pipeline_params)
    #
    # # End of map
    # end = time.time()
    # map_time = end - start
    #
    # # Delete intermediate files
    # af.delete_files(storage, pipeline_params,
    #                 cloud_prefixes=[f'map_index_files/{pipeline_params.execution_name}/',
    #                                 f'corrected_index/{pipeline_params.execution_name}/',
    #                                 f'filtered_map_files/{pipeline_params.execution_name}/'],
    #                 local_files=[f'/tmp/{pipeline_params.execution_name}/{map1_cachefile}',
    #                              f'/tmp/{pipeline_params.execution_name}/{map2_cachefile}',
    #                              f'/tmp/{pipeline_params.execution_name}/{correction_cachefile}'])
    #
    # return map_time
