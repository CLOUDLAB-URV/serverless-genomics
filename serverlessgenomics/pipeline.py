import shelve
import logging

import lithops

from .mapping.map_caller import run_full_alignment
from .preprocessing.preprocess_fasta import prepare_fasta_chunks
from .preprocessing.preprocess_fastq import prepare_fastq_chunks

# map/reduce functions and executor
from .cachedlithops import CachedLithopsInvoker

from .parameters import PipelineRun, Lithops, validate_parameters
from .constants import CACHE_PATH
from .utils import setup_logging, log_parameters

logger = logging.getLogger(__name__)


class VariantCallingPipeline:
    def __init__(self, override_id=None, **params):
        params['override_id'] = override_id
        self.parameters: PipelineRun = validate_parameters(params)
        self.fastq_chunks = None
        self.fasta_chunks = None
        self.alignment_batches = None
        self._setup()

    def _setup(self):
        setup_logging(self.parameters.log_level)
        logger.info('Init Serverless Variant Calling Pipeline')

        if self.parameters.log_level == logging.DEBUG:
            log_parameters(self.parameters)

        self.lithops = Lithops(storage=lithops.storage.Storage(), invoker=CachedLithopsInvoker(self.parameters))

        with shelve.open(CACHE_PATH) as cache:
            cache[f'{self.parameters.run_id}/parameters'] = self.parameters

    @classmethod
    def restore_run(cls, run_id: str):
        self = cls.__new__(cls)
        key = f'{run_id}/parameters'
        with shelve.open(CACHE_PATH) as cache:
            if key not in cache:
                raise KeyError(f'run {run_id} not found in cache')
            self.parameters = cache[key]
        self._setup()
        return self

    def preprocess(self):
        """
        Prepare requested input data for alignment
        """
        self.fastq_chunks = prepare_fastq_chunks(self.parameters, self.lithops)
        # fetch_fastq_chunk(self.fastq_chunks[0], 'test.fastq', self.lithops.storage, self.parameters.fastq_path, self.parameters.storage_bucket, self.parameters.fastqgz_idx_keys[0])
        self.fasta_chunks = prepare_fasta_chunks(self.parameters, self.lithops)
        # fetch_fasta_chunk(self.fasta_chunks[0], 'test', self.lithops.storage, self.parameters.fasta_path)

    def align_reads(self):
        """
        Alignment map pipeline step
        """
        assert self.fasta_chunks is not None and self.fastq_chunks is not None, 'generate chunks first!'
        run_full_alignment(self.parameters, self.lithops, self.fasta_chunks, self.fastq_chunks)

    # TODO implement reduce stage
    def reduce(self):
        raise NotImplementedError()

    def run_pipeline(self):
        """
        Execute all pipeline steps in order
        """
        self.preprocess()
        self.align_reads()

    def clean_all(self):
        keys = self.lithops.storage.list_keys(self.parameters.storage_bucket, prefix=self.parameters.fastqgz_idx_prefix)
        self.lithops.storage.delete_objects(self.parameters.storage_bucket, keys)

        keys = self.lithops.storage.list_keys(self.parameters.storage_bucket, prefix=self.parameters.faidx_prefix)
        self.lithops.storage.delete_objects(self.parameters.storage_bucket, keys)

        keys = self.lithops.storage.list_keys(self.parameters.storage_bucket, prefix=self.parameters.tmp_prefix)
        self.lithops.storage.delete_objects(self.parameters.storage_bucket, keys)
