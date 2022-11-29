import time
import sys
from random import randint
from .preprocessing import prepare_fasta, prepare_fastq, GenerateAlignmentIterdata
from .parameters import PipelineParameters
import pathlib
import logging
from lithops import Storage

# map/reduce functions and executor
from . import map_reduce_caller
from .alignment_mapper import AlignmentMapper
from . import metadata as sra_meta

from .parameters import validate_parameters, PipelineParameters
from .utils import setup_logging, log_parameters

logger = logging.getLogger(__name__)


class VariantCallingPipeline:
    def __init__(self, **kwargs):
        self.parameters: PipelineParameters = validate_parameters(kwargs)
        setup_logging(self.parameters.log_level)
        logger.info('Init Serverless Variant Calling Pipeline')
        if self.parameters.log_level == logging.DEBUG:
            log_parameters(self.parameters)

    def preprocess(self):
        storage = Storage()
        size_chunk_w = int(
            int(storage.head_object(self.parameters.bucket, 
                                    self.parameters.fasta_folder + self.parameters.fasta_file)['content-length']) / int(self.parameters.fasta_workers))

        ###################################################################
        #### GENERATE LIST OF FASTQ CHUNKS (BYTE RANGES)
        ###################################################################

        num_spots = 0
        metadata = sra_meta.SraMetadata()
        arr_seqs = [self.parameters.fq_seqname]
        if self.parameters.datasource == "SRA":
            accession = metadata.efetch_sra_from_accessions(arr_seqs)
            seq_type = accession['pairing'].to_string(index=False)
            num_spots = accession['spots'].to_string(index=False)
            fastq_size = accession['run_size'].to_string(index=False)
            print("Retrieving data from sra...")
            print("Sequence type: " + seq_type)
            print("Number of spots: " + num_spots)
            print("fastq size: " + fastq_size)

        fastq_list = prepare_fastq(self.parameters.fastq_read_n, self.parameters.fq_seqname, num_spots)

        ###################################################################
        #### GENERATE LIST OF FASTA CHUNKS (if not present)
        ###################################################################

        fasta_index = prepare_fasta(self.parameters)

        ###################################################################
        #### GENERATE ITERDATA AND PREPROCESSING SUMMARY
        ###################################################################
        # The iterdata consists of an array where each element is a pair of a fastq chunk and a fasta chunk.
        # Since each fastq chunk needs to be paired with each fasta chunk, the total number of elements in
        # iterdata will be n_fastq_chunks * n_fasta_chunks.

        generate_aligment = GenerateAlignmentIterdata(self.parameters, fastq_list, fasta_index,
                                                      self.parameters.fasta_folder + self.parameters.fasta_file)
        return generate_aligment()

    # TODO implement alignment stage
    def map_alignment(self):
        ...

    # TODO implement reduce stage
    def reduce(self):
        ...

    def run_pipeline(self):
        """
        Execute all pipeline steps in order
        """
        iterdata, num_chunks = self.preprocess()
        

        if not self.parameters.pre_processing_only:
            ###################################################################
            #### MAP-REDUCE
            ###################################################################
            mapfunc = AlignmentMapper(pathlib.Path(self.parameters.fasta_file).stem, self.parameters)
            map_time = map_reduce_caller.map_reduce(self.parameters, iterdata, mapfunc, num_chunks)

            ###################################################################
            #### MAP-REDUCE SUMMARY
            ###################################################################
            print("MAP-REDUCE SUMMARY")
            print("map phase: execution_time_total_varcall: " + str(map_time) + "s")
