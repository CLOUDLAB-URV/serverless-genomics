import time
import re
import sys
from random import randint
from . import fastq_functions as fq_func
from . import fasta_functions as fa_func
from .fastaPartitionerIndex import FunctionsFastaIndex
from .parameters import PipelineParameters
import pathlib
import logging
from lithops import Storage
from fastaPartitionerIndex import FunctionsFastaIndex

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

    # TODO move this method to preprocessing module
    def generate_alignment_iterdata(self, args: PipelineParameters, list_fastq: list, fasta_index: str, fasta_file_path,
                                    iterdata_n):
        """
        Creates the lithops iterdata from the fasta and fastq chunk lists
        """
        # Number of fastq chunks processed. If doing a partial execution, iterdata_n will need to be multiple of the number of fasta chunks
        # so the index correction is done properly.
        num_chunks = 0
        # TODO this try/except block is reduntant - an exception already exits the program
        try:
            functions = FunctionsFastaIndex(fasta_index, fasta_file_path)
            fasta_chunks = functions.get_chunks(args)

            iterdata = []
            for fastq_key in list_fastq:
                num_chunks += 1
                for i, chunk in enumerate(fasta_chunks):
                    iterdata.append({'fasta_chunk': {'key_fasta': fasta_file_path, 'key_index': fasta_index, 'id': i,
                                                     'chunk': chunk}, 'fastq_chunk': fastq_key,
                                     'exec_param': args.execution_name})

                    # Limit the length of iterdata if iterdata_n is not null.
            if iterdata_n is not None and iterdata_n != "all":
                iterdata = iterdata[0:int(iterdata_n)]
                if len(iterdata) % len(fasta_chunks) != 0:
                    raise Exception(
                        f"ERROR. Number of elements in iterdata must be multiple of the number of fasta chunks (iterdata: {len(iterdata)}, data_index: {len(fasta_chunks)}).")
                else:
                    num_chunks = int(
                        re.sub('^[\s|\S]*number\':\s(\d*),[\s|\S]*$', r"\1", str(iterdata[-1]['fastq_chunk'])))
            if not iterdata:
                raise Exception('ERROR. Iterdata not generated')
        except Exception as e:
            print(e)
            exit()

        return iterdata, len(fasta_chunks), num_chunks

    # TODO implement preprocess stage
    def preprocess(self):
        ...

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
        # TODO move to preprocessing stage
        storage = Storage()
        size_chunk_w = int(
            int(storage.head_object(self.parameters.bucket, self.parameters.fasta_folder + self.parameters.fasta_file)[
                    'content-length']) / int(
                self.parameters.fasta_workers))

        run_id = str(randint(1000, 9999))

        pp_start = time.time()

        # Run settings summary
        # print("Variant Caller - Cloudbutton Genomics Use Case demo")
        # print("starting pipeline at " + str(time.time()) + " - " + str(
        #     time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())))
        #
        # print("run id: " + run_id)
        #
        # print("RUN SETTINGS")
        # print("command line: ")
        # print(sys.argv[0] + str(sys.argv))
        # print("\nBucket name: %s\n" % args.bucket)
        # print("Sequencing type: " + args.seq_type + "\n")
        #
        # print("INPUT FILES")
        # print("Fastq file 1: %s" % args.fastq_file)
        # print("Fastq file 2: %s" % args.fastq_file2)
        # print("Fasta file: %s" % args.fasta_file)
        #
        # print("\nFILE SPLITTING SETTINGS")
        # print("Fastq chunk size: %s lines" % str(args.fastq_chunk_size))
        # print("Fasta number of workers: %s" % str(args.fasta_workers))
        # print("Chunk per workers: %s B" % str(size_chunk_w))
        #
        # print("\nOTHER RUN SETTINGS")
        # if args.function_n == "all":
        #     print("Number of functions spawned: %s" % args.function_n)
        # else:
        #     print("Number of functions spawned: %s" % args.iterdata_n)
        # print("Runtime used: %s" % args.runtime_id)
        # print("Runtime memory - map function: %s" % args.runtime_mem)
        # print("Runtime memory - reduce function: %s" % args.runtime_mem_r)
        # print("Runtime storage - map function: %s" % args.runtime_storage)
        # print("Reduce load balancer method: %s" % args.lb_method)
        # print("############\n")

        ###################################################################
        #### 1. GENERATE LIST OF FASTQ CHUNKS (BYTE RANGES)
        ###################################################################
        t0 = time.time()

        # TODO move to preprocessing stage
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

        # TODO move to preprocessing stage
        # Generate fastq chunks
        fastq_list = fq_func.prepare_fastq(self.parameters.fastq_read_n, self.parameters.fq_seqname, num_spots)

        t1 = time.time()
        print(f'PP:0: fastq list: execution_time: {t1 - t0}: s')

        ###################################################################
        #### 2. GENERATE LIST OF FASTA CHUNKS (if not present)
        ###################################################################
        t2 = time.time()

        # Fasta File Chunk Prefix:
        # the prefix contains the name of the fasta reference minus the suffix,
        # followed by block length and then "_split_" i.e. for hg19.fa
        # with chunks with block length of 400000 it would be hg19_400000_split_
        # TODO move to preprocessing stage
        fasta_chunks_prefix = pathlib.Path(self.parameters.fasta_file).stem

        # Generate fasta chunks
        fasta_index = fa_func.prepare_fasta(self.parameters)

        t3 = time.time()
        print(f'PP:0: fasta list: execution_time: {t3 - t2}: s')

        ###################################################################
        #### 3. GENERATE ITERDATA
        ###################################################################
        # The iterdata consists of an array where each element is a pair of a fastq chunk and a fasta chunk.
        # Since each fastq chunk needs to be paired with each fasta chunk, the total number of elements in
        # iterdata will be n_fastq_chunks * n_fasta_chunks.
        print("\nGenerating iterdata")
        # TODO move to preprocessing stage
        # Generate iterdata
        iterdata, fastq_set_n, num_chunks = self.generate_alignment_iterdata(self.parameters, fastq_list, fasta_index,
                                                                             self.parameters.fasta_folder + self.parameters.fasta_file,
                                                                             self.parameters.iterdata_n)

        iterdata_n = len(iterdata)
        ###################################################################
        #### 4. PREPROCESSING SUMMARY
        ###################################################################
        print("\nITERDATA LIST")
        print("number of fasta chunks: " + str(fastq_set_n))
        print("number of fastq chunks: " + str(len(fastq_list)))
        print("fasta x fastq chunks: " + str(fastq_set_n * len(fastq_list)))

        print("number of iterdata elements: " + str(iterdata_n))
        pp_end = time.time()
        print(f'PP:0: preprocessing: execution_time: {pp_end - pp_start}: s')

        if not self.parameters.pre_processing_only:
            ###################################################################
            #### 5. MAP-REDUCE
            ###################################################################
            mapfunc = AlignmentMapper(fasta_chunks_prefix, self.parameters)
            map_time = map_reduce_caller.map_reduce(self.parameters, iterdata, mapfunc, num_chunks)

            ###################################################################
            #### 6. MAP-REDUCE SUMMARY
            ###################################################################
            print("MAP-REDUCE SUMMARY")
            print("map phase: execution_time_total_varcall: " + str(map_time) + "s")
