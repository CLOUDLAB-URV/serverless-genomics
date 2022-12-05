import time
import sys
from random import randint
from .preprocessing import prepare_fasta, prepare_fastq, generate_alignment_iterdata
from .parameters import PipelineParameters
import pathlib
import logging

# map/reduce functions and executor
from .mapping import map_caller
from .mapping.alignment_mapper import AlignmentMapper
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
        ###################################################################
        #### GENERATE LIST OF FASTQ CHUNKS (BYTE RANGES)
        ###################################################################

        if self.parameters.datasource == "SRA":
            metadata = sra_meta.SraMetadata()
            accession = metadata.efetch_sra_from_accessions([self.parameters.fq_seqname])
            seq_type = accession['pairing'].to_string(index=False)
            num_spots = accession['spots'].to_string(index=False)
            fastq_size = accession['run_size'].to_string(index=False)
            print("Retrieving data from sra...")
            print("Sequence type: " + seq_type)
            print("Number of spots: " + num_spots)
            print("fastq size: " + fastq_size)
        else:
            num_spots = 0

        fastq_list = prepare_fastq(int(self.parameters.fastq_read_n), self.parameters.fq_seqname, int(num_spots))

        ###################################################################
        #### GENERATE LIST OF FASTA CHUNKS (if not present)
        ###################################################################

        # Generate the index file from the fasta file source if it does not exist and return the path in the storage to the fasta file
        fasta_index = prepare_fasta(self.parameters)

        ###################################################################
        #### GENERATE ITERDATA AND PREPROCESSING SUMMARY
        ###################################################################

        # Creates the lithops iterdata from the fasta and fastq chunk lists
        return generate_alignment_iterdata(self.parameters, fastq_list, fasta_index, self.parameters.fasta_folder + self.parameters.fasta_file)


    def map_alignment(self, iterdata, num_chunks):
        ###################################################################
        #### MAP-REDUCE
        ###################################################################
        mapfunc = AlignmentMapper(pathlib.Path(self.parameters.fasta_file).stem, self.parameters)
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
            map_time = self.map_alignment(iterdata, num_chunks)
            print("map-reduce phase execution time: " + str(map_time) + "s")