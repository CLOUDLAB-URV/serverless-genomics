import time
import sys
from random import randint

import lithops

from .preprocessing.fasta_functions import prepare_fasta
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

    def preprocess(self):
        """
        Prepare requested input data for alignment
        """
        fastq_list = prepare_fastq_chunks(self.parameters, self.lithops)
        # fasta_index = prepare_fasta(self.parameters)
        #
        # iter_data = generate_alignment_iterdata(self.parameters, fastq_list, fasta_index,
        #                                         self.parameters.fasta_folder + self.parameters.fasta_path)

    def run_alignment(self):
        ###################################################################
        #### MAP-REDUCE
        ###################################################################
        mapfunc = AlignmentMapper(pathlib.Path(self.parameters.fasta_path).stem, self.parameters)
        map_time = map_caller.map(self.parameters, iterdata, mapfunc, num_chunks)
        return map_time

    # TODO implement reduce stage
    def reduce(self):
        ...

    def run_pipeline(self):
        """
        Execute all pipeline steps in order
        """
        iterdata, num_chunks = self.preprocess()

        if not self.parameters.pre_processing_only:
            map_time = self.run_alignment(iterdata, num_chunks)
            print("map-reduce phase execution time: " + str(map_time) + "s")
