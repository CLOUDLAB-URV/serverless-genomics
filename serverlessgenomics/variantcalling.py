import logging

import lithops

from .mapping.map_caller import run_full_alignment
from .preprocessing import (
    prepare_fastq_chunks,
    prepare_fasta_chunks,
    prepare_gem_chunks,
)
from .reducer.reduce_caller import run_reducer
from .stats import Stats

from .pipeline import (
    PipelineParameters,
    PipelineRun,
    Lithops,
    validate_parameters,
    new_pipeline_run,
)
from .utils import setup_logging, get_storage_tmp_prefix
from .lithopswrapper import LithopsInvokerWrapper

logger = logging.getLogger(__name__)


class VariantCallingPipeline:
    def __init__(self, **parameters):
        run_id = parameters.pop("run_id", None)
        self.parameters: PipelineParameters = validate_parameters(parameters)
        setup_logging(self.parameters.log_level)

        logger.info("Init Serverless Variant Calling Pipeline")
        self.state: PipelineRun = new_pipeline_run(self.parameters, run_id)
        self.global_stat = Stats()

        invoker = LithopsInvokerWrapper(self.parameters.lithops_settings)
        storage = lithops.Storage()
        self.lithops = Lithops(storage=storage, invoker=invoker)

    def preprocess(self):
        """
        Prepare requested input data for alignment
        """
        with self.global_stat.timeit("prepare_fastq_chunks"):
            fastq_chunks = prepare_fastq_chunks(self.parameters, self.lithops)
        self.state.fastq_chunks = fastq_chunks

        with self.global_stat.timeit("prepare_fasta_chunks"):
            fasta_chunks = prepare_fasta_chunks(self.parameters, self.lithops)
        self.state.fasta_chunks = fasta_chunks

        with self.global_stat.timeit("prepare_gem_chunks"):
            self.state.gem_chunk_ids, gem_stats = prepare_gem_chunks(
                self.parameters, self.state.fasta_chunks, self.lithops
            )
        self.global_stat.set_value("gem_stats", [s.dump_dict() for s in gem_stats])

    def alignment(self):
        """
        Alignment map pipeline step
        """
        with self.global_stat.timeit("run_full_alignment"):
            stats = run_full_alignment(self.parameters, self.state, self.lithops)
        self.global_stat.set_value("alignment_stats", stats.dump_dict())

    def reduce(self):
        with self.global_stat.timeit("run_reducer"):
            if self.state.aligned_mpileups is None:
                mpileups_prefix = get_storage_tmp_prefix(self.state.run_id, "filtered_index_to_mpileup")
                objects = self.lithops.storage.list_objects(
                    bucket=self.parameters.storage_bucket, prefix=mpileups_prefix
                )
                keys = [obj["Key"] for obj in objects]
                self.state.aligned_mpileups = {i: key for i, key in enumerate(keys)}
            stats = run_reducer(
                self.parameters,
                self.state,
                self.lithops,
            )
        self.global_stat.set_value("reduce_stats", stats.dump_dict())

    def pipeline_stats(self):
        stats, params = PipelineRunStats(), PipelineRunStats()

        stats.store_size_data("fasta_path", str(self.parameters.fasta_path))
        stats.store_size_data("fastq_path", str(self.parameters.fastq_path))
        stats.store_size_data("fastq_chunks", self.parameters.fastq_chunks)
        stats.store_size_data("fasta_chunks", self.parameters.fasta_chunks)
        stats.store_size_data("run_id", str(self.parameters.run_id))
        if self.parameters.fastq_chunk_range is not None:
            stats.store_size_data("fastq_range", str(self.parameters.fastq_chunk_range))

        stats.store_dictio(params.get_stats(), "pipeline_params")
        return stats

    def run_pipeline(self):
        """
        Execute all pipeline steps in order
        """
        self.preprocess()
        self.alignment()
        self.reduce()

    def clean_temp_data(self):
        prefix = get_storage_tmp_prefix(self.state.run_id, "")
        logger.info("Going to delete all temporary data with prefix %s", prefix)
        keys = self.lithops.storage.list_keys(self.parameters.storage_bucket, prefix=prefix)
        self.lithops.storage.delete_objects(self.parameters.storage_bucket, keys)

    def clean_all(self):
        logger.info("Going to delete all FASTQGZ Indexes")
        keys = self.lithops.storage.list_keys(
            self.parameters.storage_bucket, prefix=self.parameters.fastqgz_idx_prefix
        )
        self.lithops.storage.delete_objects(self.parameters.storage_bucket, keys)

        logger.info("Going to delete all FAIDX Indexes")
        keys = self.lithops.storage.list_keys(self.parameters.storage_bucket, prefix=self.parameters.faidx_prefix)
        self.lithops.storage.delete_objects(self.parameters.storage_bucket, keys)

        logger.info("Going to delete all GEM Indexes")
        keys = self.lithops.storage.list_keys(self.parameters.storage_bucket, prefix=self.parameters.gem_index_prefix)
        self.lithops.storage.delete_objects(self.parameters.storage_bucket, keys)
