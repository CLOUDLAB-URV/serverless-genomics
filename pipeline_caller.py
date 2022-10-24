import argparse
import os.path
import time
import re
import sys
from random import randint
import Metadata
import fastq_functions as fq_func
import fasta_functions as fa_func
from varcall_arguments import Arguments

# map/reduce functions and executor
from map_reduce_executor import MapReduce
from map_functions import MapFunctions

class PipelineCaller:
    def generate_alignment_iterdata(self, list_fastq, list_fasta, iterdata_n):
        """
        Creates the lithops iterdata from the fasta and fastq chunk lists
        """
        iterdata = []
        
        # Number of fastq chunks processed. If doing a partial execution, iterdata_n will need to be multiple of the number of fasta chunks
        # so the index correction is done properly.
        num_chunks = 0  
        
        # Generate iterdata
        for fastq_key in list_fastq:
            num_chunks += 1
            for fasta_key in list_fasta:
                iterdata.append({'fasta_chunk': fasta_key, 'fastq_chunk': fastq_key})

        # Limit the length of iterdata if iterdata_n is not null.
        if iterdata_n is not None or iterdata_n == "all": 
            iterdata = iterdata[0:int(iterdata_n)]
            if(len(iterdata)%len(list_fasta)!=0):
                raise Exception("If ")
            else: 
                num_chunks = len(iterdata)//len(list_fasta)

        return iterdata, num_chunks
    
    def __call__(self, args: Arguments):
        ###################################################################
        #### START THE PIPELINE
        ###################################################################
        run_id=str(randint(1000,9999))
        stage="PP"  
        id="X"      # this is for function id, which is not present in this case
        
        PP_start = time.time()

        # Run settings summary
        print("Variant Caller - Cloudbutton Genomics Use Case demo")
        print("starting pipeline at "+ str(time.time()) + " - " + str(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())))
        
        print("run id: "+run_id)
        
        print("RUN SETTINGS")
        print("command line: ")
        print(sys.argv[0]+str(sys.argv))
        print("\nBucket name: %s\n" % args.bucket )
        print("Sequencing type: " + args.seq_type + "\n")

        print("INPUT FILES")
        print("Fastq file 1: %s" % args.fastq_file )
        print("Fastq file 2: %s" % args.fastq_file2 )
        print("Fasta file: %s" % args.fasta_file )

        print("\nFILE SPLITTING SETTINGS")
        print("Fastq chunk size: %s lines" % str(args.fastq_chunk_size) )
        print("Fasta chunk size: %s characters" % str(args.fasta_chunk_size) )
        print("Fasta overlap size: %s characters" % str(args.fasta_char_overlap) )

        print("\nOTHER RUN SETTINGS")
        if args.function_n == "all":
            print("Number of functions spawned: %s" % args.function_n )
        else:
            print("Number of functions spawned: %s" % args.iterdata_n )
        print("Runtime used: %s" % args.runtime_id )
        print("Runtime memory - map function: %s" % args.runtime_mem )
        print("Runtime memory - reduce function: %s" % args.runtime_mem_r )
        print("Runtime storage - map function: %s" % args.runtime_storage )
        print("Reduce load balancer method: %s" % args.lb_method )
        print("############\n")


        ###################################################################
        #### 1. GENERATE LIST OF FASTQ CHUNKS (BYTE RANGES)
        ###################################################################
        t0 = time.time()

        num_spots = 0
        metadata = Metadata.SraMetadata()
        arr_seqs = [args.fq_seqname]
        if args.datasource == "SRA":
            accession = metadata.efetch_sra_from_accessions(arr_seqs)
            seq_type = accession['pairing'].to_string(index=False)
            num_spots = accession['spots'].to_string(index=False)
            fastq_size = accession['run_size'].to_string(index=False)
            print("Retrieving data from sra...")
            print("Sequence type: " + seq_type)
            print("Number of spots: " + num_spots)
            print("fastq size: " + fastq_size)

        # Generate fastq chunks
        fastq_list = fq_func.prepare_fastq(args.fastq_read_n, args.fq_seqname, num_spots)

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
        fasta_chunks_prefix = re.sub("\.fasta|\.fa|\.fas", "_" + str(args.fasta_chunk_size) + "split_", args.fasta_file)
        
        # Generate fasta chunks
        fasta_list = fa_func.prepare_fasta(args.cloud_adr, args.runtime_id, args.fasta_bucket, args.fasta_folder, args.fasta_file, args.split_fasta_folder, args.fasta_chunk_size, fasta_chunks_prefix, args.fasta_char_overlap)
        
        t3 = time.time()
        print(f'PP:0: fasta list: execution_time: {t3 - t2}: s')

        
        ###################################################################
        #### 3. GENERATE ITERDATA
        ###################################################################
        # The iterdata consists of an array where each element is a pair of a fastq chunk and a fasta chunk.
        # Since each fastq chunk needs to be paired with each fasta chunk, the total number of elements in
        # iterdata will be n_fastq_chunks * n_fasta_chunks.
        
        print("\nGenerating iterdata")
        t4 = time.time()
        
        # Generate iterdata
        iterdata, num_chunks = self.generate_alignment_iterdata(fastq_list, fasta_list, args.iterdata_n)
        
        iterdata_n=len(iterdata)
        fastq_set_n=len(fasta_list) # number of fastq files aligned to each fasta chunk.
        
        
        ###################################################################
        #### 4. PREPROCESSING SUMMARY
        ###################################################################
        print("\nITERDATA LIST")
        print("number of fasta chunks: " + str(fastq_set_n))
        print("number of fastq chunks: " + str(len(fastq_list)))
        print("fasta x fastq chunks: "+ str(len(fasta_list)*len(fastq_list)))
        print("number of iterdata elements: " + str(iterdata_n))
        PP_end = time.time()
        print(f'PP:0: preprocessing: execution_time: {PP_end - PP_start}: s')


        if args.pre_processing_only==False:
            ###################################################################
            #### 5. MAP-REDUCE
            ###################################################################
            stage="MR"
            id=0
            
            mapfunc = MapFunctions(args.fq_seqname, args.datasource, seq_type, args.debug, args.bucket, args.fastq_folder, args.idx_folder, fasta_chunks_prefix, args.fasta_bucket, stage, args.file_format, args.tolerance)

            mapreduce = MapReduce(mapfunc, args.runtime_id, args.runtime_mem, args.runtime_mem_r, args.buffer_size, args.file_format, args.func_timeout_map, args.func_timeout_reduce, 'DEBUG', args.bucket, stage, id, args.debug, args.skip_map, args.lb_method, iterdata_n, num_chunks, args.fq_seqname)
            map_time = mapreduce(iterdata)

            
            ###################################################################
            #### 6. MAP-REDUCE SUMMARY
            ###################################################################
            print("MAP-REDUCE SUMMARY")
            print("map phase: execution_time_total_varcall: "+str(map_time)+"s")