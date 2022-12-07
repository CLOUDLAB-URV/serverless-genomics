import time
import sys
from random import randint

import lithops

from .preprocessing.preprocess_fasta import prepare_fasta_chunks
from .preprocessing.alignment_iterdata import generate_alignment_iterdata
from .preprocessing.preprocess_fastq import prepare_fastq_chunks
import pathlib
import logging

# map/reduce functions and executor
from .mapping import map_caller
from .mapping.alignment_mapper import AlignmentMapper
from .lithopsproxy import LithopsProxy

from .parameters import PipelineRun, Lithops, validate_parameters
from .utils import setup_logging, log_parameters, S3Path

logger = logging.getLogger(__name__)


class VariantCallingPipeline:
    def __init__(self, **params):
        self.parameters: PipelineRun = validate_parameters(params)

        setup_logging(self.parameters.log_level)
        logger.info('Init Serverless Variant Calling Pipeline')

        if self.parameters.log_level == logging.DEBUG:
            log_parameters(self.parameters)

        self.lithops = Lithops(storage=lithops.storage.Storage(), invoker=LithopsProxy())
        self.fastq_chunks = None
        self.fasta_chunks = None

    def preprocess(self):
        """
        Prepare requested input data for alignment
        """
        self.fastq_chunks = prepare_fastq_chunks(self.parameters, self.lithops)
        self.fasta_chunks = prepare_fasta_chunks(self.parameters, self.lithops)

    def align_reads(self):
        """
        Alignment map pipeline step
        """
        assert self.fasta_chunks is not None and self.fastq_chunks is not None, 'generate chunks first!'

        mapfunc = AlignmentMapper(pathlib.Path(self.parameters.fasta_path).stem, self.parameters)
        map_time = map_caller.map(self.parameters, iterdata, mapfunc, num_chunks)
        return map_time

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