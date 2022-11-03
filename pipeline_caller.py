import time
import re
import sys
from random import randint
import metadata as Metadata
import fastq_functions as fq_func
import fasta_functions as fa_func
from varcall_arguments import Arguments
import pathlib
from lithops import Storage

# map/reduce functions and executor
import map_reduce_caller as map_reduce
from alignment_mapper import AlignmentMapper

class PipelineCaller:
    def __generate_info_seq(self, args, storage, fasta_file_path, values, data_index, i):
        try: 
            next_val = data_index[i+1].split(' ')
            last_byte_plus = int(next_val[1])-1 if next_val[0] != values[0] else int(next_val[2])-1
        except:
            last_byte_plus = int(storage.head_object(args.bucket, fasta_file_path)['content-length']) - 1
        return {'offset_head': int(values[1]), 'offset_base':  int(values[2]), 'last_byte+': last_byte_plus}

    def generate_alignment_iterdata(self, args: Arguments, list_fastq: list, fasta_index: str, fasta_file_path, iterdata_n):
        """
        Creates the lithops iterdata from the fasta and fastq chunk lists
        """
        storage = Storage()  
        
        # Number of fastq chunks processed. If doing a partial execution, iterdata_n will need to be multiple of the number of fasta chunks
        # so the index correction is done properly.
        num_chunks = 0  
        data_index = []
        iterdata = []
        fasta_chunks = []
        try:
            fasta = storage.head_object(args.bucket, args.fasta_folder+args.fasta_file)
            fa_chunk = int(int(fasta['content-length']) / int(args.fasta_workers))
            data_index = storage.get_object(args.bucket, fasta_index).decode('utf-8').split('\n')
                          
            values = data_index[0].split(' ')
            last_seq = values[4]
            tmp_seq = self.__generate_info_seq(args, storage, fasta_file_path, values, data_index, 0)
            fa_chunk = [tmp_seq] 
            for i, sequence in enumerate(data_index[1:], start=1):                    
                values = sequence.split(' ')
                if last_seq == values[4]:
                    tmp_seq = self.__generate_info_seq(args, storage, fasta_file_path, values, data_index, i)
                else:
                    fa_chunk.append(tmp_seq)
                    fasta_chunks.append(fa_chunk)
                    tmp_seq = self.__generate_info_seq(args, storage, fasta_file_path, values, data_index, i)
                    fa_chunk = [tmp_seq]
                last_seq = values[4]
                
            for fastq_key in list_fastq:
                num_chunks += 1
                for i, chunk in enumerate(fasta_chunks):
                    iterdata.append({'fasta_chunk': {'key_fasta': fasta_file_path, 'key_index': fasta_index, 'id': i, 'chunk': chunk}, 'fastq_chunk': fastq_key}) 

            # Limit the length of iterdata if iterdata_n is not null.
            if iterdata_n is not None and iterdata_n != "all":
                iterdata = iterdata[0:int(iterdata_n)]
                if(len(iterdata)%len(fasta_chunks)!=0):
                    raise Exception(f"ERROR. Number of elements in iterdata must be multiple of the number of fasta chunks (iterdata: {len(iterdata)}, data_index: {len(fasta_chunks)}).")
                else:
                    num_chunks = int(re.sub('^[\s|\S]*number\':\s(\d*),[\s|\S]*$', r"\1", str(iterdata[-1]['fastq_chunk'])))
            if not iterdata:
                raise Exception('ERROR. Iterdata not generated')
        except Exception as e:
            print(e)
            exit()

        return iterdata, len(fasta_chunks), num_chunks
    
    def execute_pipeline(self, args: Arguments):
        ###################################################################
        #### START THE PIPELINE
        ###################################################################
        storage = Storage()
        size_chunk_w = int(int(storage.head_object(args.bucket, args.fasta_folder+args.fasta_file)['content-length']) / int(args.fasta_workers))

        run_id=str(randint(1000,9999))
        
        pp_start = time.time()

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
        print("Fasta number of workers: %s" % str(args.fasta_workers) ) 
        print("Chunk per workers: %s B" % str(size_chunk_w) ) 
        
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
        fasta_chunks_prefix = pathlib.Path(args.fasta_file).stem
        
        # Generate fasta chunks
        fasta_index = fa_func.prepare_fasta(args)
        
        t3 = time.time()
        print(f'PP:0: fasta list: execution_time: {t3 - t2}: s')

        
        ###################################################################
        #### 3. GENERATE ITERDATA
        ###################################################################
        # The iterdata consists of an array where each element is a pair of a fastq chunk and a fasta chunk.
        # Since each fastq chunk needs to be paired with each fasta chunk, the total number of elements in
        # iterdata will be n_fastq_chunks * n_fasta_chunks.
        print("\nGenerating iterdata")
        
        # Generate iterdata
        iterdata, fastq_set_n, num_chunks = self.generate_alignment_iterdata(args, fastq_list, fasta_index, args.fasta_folder+args.fasta_file, args.iterdata_n)

        iterdata_n=len(iterdata)
        ###################################################################
        #### 4. PREPROCESSING SUMMARY
        ###################################################################
        print("\nITERDATA LIST")
        print("number of fasta chunks: " + str(fastq_set_n))
        print("number of fastq chunks: " + str(len(fastq_list)))
        print("fasta x fastq chunks: "+ str(fastq_set_n*len(fastq_list)))

        print("number of iterdata elements: " + str(iterdata_n))
        pp_end = time.time()
        print(f'PP:0: preprocessing: execution_time: {pp_end - pp_start}: s')


        if args.pre_processing_only==False:
            ###################################################################
            #### 5. MAP-REDUCE
            ###################################################################
            mapfunc = AlignmentMapper(fasta_chunks_prefix, args)
            map_time = map_reduce.map_reduce(args, iterdata, mapfunc, num_chunks)
            
            ###################################################################
            #### 6. MAP-REDUCE SUMMARY
            ###################################################################
            print("MAP-REDUCE SUMMARY")
            print("map phase: execution_time_total_varcall: "+str(map_time)+"s")
