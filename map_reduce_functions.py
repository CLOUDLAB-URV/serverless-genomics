# PACKAGES
# generic packages
import os.path
import re
from random import randint

# specific packages (cloud computing, parallelisation etc)
import multiprocessing
import subprocess as sp
import pandas as pd
from numpy import int64


# demo functions and packages
import lithopsgenetics
import lithopsgenetics.auxiliaryfunctions as af
import math
import boto3
import time


class MapReduceFunctions:
    def __init__(self, fq_seqname, datasource, seq_type, debug, BUCKET_NAME, fastq_folder, idx_folder,
                  fasta_chunks_prefix, FASTA_BUCKET, gem_test, stage, file_format, tolerance):
        self.fq_seqname = fq_seqname
        self.datasource = datasource
        self.seq_type = seq_type
        self.debug = debug
        self.BUCKET_NAME = BUCKET_NAME
        self.fastq_folder = fastq_folder
        self.idx_folder = idx_folder
        self.fasta_chunks_prefix = fasta_chunks_prefix
        self.FASTA_BUCKET = FASTA_BUCKET
        self.gem_test = gem_test
        self.stage = stage 
        self.file_format = file_format
        self.tolerance = tolerance
    
    ###################################################################
    ###################################################################
    #### MAP FUNCTIONS
    ###################################################################
    ###################################################################

    def map_alignment1(self, id, fasta_chunk, fastq_chunk, storage):
        """
        First map function to filter and map fasta + fastq chunks. Some intermediate files
        are uploaded into the cloud storage for subsequent index correction, after which
        the final part of the map function (map_alignment2) can be executed
        """
        
        #Global variables imported from main file
        fq_seqname = self.fq_seqname
        datasource = self.datasource
        seq_type = self.seq_type
        debug = self.debug
        BUCKET_NAME = self.BUCKET_NAME
        fastq_folder = self.fastq_folder
        idx_folder = self.idx_folder
        fasta_chunks_prefix = self.fasta_chunks_prefix
        FASTA_BUCKET = self.FASTA_BUCKET
        gem_test = self.gem_test
        
        # Settings summary    
        stage="A"+str(randint(1000,9999))
        cpus=multiprocessing.cpu_count()
        
        # CONTROL VARIABLES 
        # fasta chunk ID
        fasta_n=re.sub(r'^\S*split_0*(\d*)\S*fasta', r'\1', fasta_chunk)
        fastq_n=re.sub('^[\s|\S]*number\':\s(\d*),[\s|\S]*$', r"\1", str(fastq_chunk))
        func_key=fq_seqname+"_fq"+fastq_n+"-fa"+fasta_n
        
        ###################################################################
        #### PROCESSING FASTQ CHUNKS
        ###################################################################
        # Download fastq chunk depending on source
        fastq1 =""
        fastq2 =""
        if datasource == "s3":
            if seq_type == "paired-end":
                af.printl("processing paired-end fastq chunks",stage,id,debug)
                fastq1 = lithopsgenetics.fastq_to_mapfun("fastq1", fastq_chunk[0][0], fastq_chunk[0][1], BUCKET_NAME, fastq_folder, idx_folder, datasource, stage,id,debug)
                fastq2 = lithopsgenetics.fastq_to_mapfun("fastq2", fastq_chunk[1][0], fastq_chunk[1][1], BUCKET_NAME, fastq_folder, idx_folder, datasource, stage,id,debug)
                base_name = os.path.splitext(fastq1)[0]#+'.pe'
            else:   # single-end sequencing
                af.printl("processing single-end fastq chunk",stage,id,debug)
                fastq1 = lithopsgenetics.fastq_to_mapfun("fastq", fastq_chunk[0], fastq_chunk[1], BUCKET_NAME, fastq_folder, idx_folder, datasource, stage,id,debug)
                base_name = os.path.splitext(fastq1)[0]+'.se'
                fastq2 = "no"      
        elif datasource == "SRA":
            if seq_type == "paired-end":
                af.printl("processing paired-end fastq chunks",stage,id,debug)
                fastq1 = lithopsgenetics.fastq_to_mapfun("fastq", fastq_chunk[0], fastq_chunk[1], BUCKET_NAME, fastq_folder, idx_folder, datasource, stage,id,debug)
                base_name = os.path.splitext(fastq1)[0]#+'.pe'
                fastq2 = "yes"
            else:   # single-end sequencing
                af.printl("processing single-end fastq chunk",stage,id,debug)
                fastq1 = lithopsgenetics.fastq_to_mapfun("fastq", fastq_chunk[0], fastq_chunk[1], BUCKET_NAME, fastq_folder, idx_folder, datasource, stage,id,debug)
                fastq2 = "no"
                base_name = os.path.splitext(fastq1)[0]#+'.se'
        
        
        ###################################################################
        #### PROCESSING FASTA CHUNKS
        ###################################################################
        # Checking if gem files exist in storage from previous runs
        gem_list = []
        gem_folder = "gem-chunks/"
        gem_ref=""
        gem_index_present=""
        try:
            gem_list = storage.list_keys(BUCKET_NAME, prefix=gem_folder + fasta_chunks_prefix)  # Gem files found
        except:
            af.printl("split gem folder empty / not found", stage, id, debug)                   # Gem files not found

        # Gem files have not been generated
        if not gem_list or gem_list == []:  
            # 2. INDEXING FASTA CHUNK
            fasta_chunk_folder_file = fasta_chunk.split("/")
            fasta = af.copy_to_runtime(storage, FASTA_BUCKET, fasta_chunk_folder_file[0]+"/", fasta_chunk_folder_file[1], stage, id, debug)
            gem_ref_nosuffix = os.path.splitext(fasta)[0]
            gem_ref = gem_ref_nosuffix + '.gem'
            
            indexer_output = sp.run(['gem-indexer', '--input', fasta, '--threads', str(cpus), '-o', gem_ref_nosuffix], capture_output=True)
        else:
            # Download gem files generated from previous runs
            gem_index_present=True
            for gem_file in gem_list:
                query=r'' + str(fasta_n) + '.gem$'
                if re.search(query, gem_file):
                    gem_file_name=fasta_chunks_prefix+str(fasta_n)+".gem"
                    gem_ref = af.copy_to_runtime(storage, BUCKET_NAME, gem_folder, gem_file_name, stage, id, debug)
                else:
                    af.printl("unmatched gem file "+gem_file)


        ###################################################################
        #### GENERATE ALIGNMENT AND ALIGNMENT INDEX (FASTQ TO MAP)
        ###################################################################
        mapper_output = sp.run(['/function/bin/map_index_and_filter_map_file_cmd_awsruntime.sh', gem_ref, fastq1, fastq2, base_name, datasource, seq_type], capture_output=True)
        
        if gem_test == True:
            af.printl("running only gem indexer and mapper - skipping rest of map function", stage, id, debug)
            return 0

        # Initialize file names
        map_index_file = base_name + "_map.index.txt"
        os.rename(map_index_file,"/tmp/" + func_key+ "_map.index.txt")
        map_index_file = "/tmp/" + func_key+ "_map.index.txt"
        old_filtered_map_file = base_name + "_filt_wline_no.map"
        filtered_map_file = base_name + "_" + str(id) + "_filt_wline_no.map"
        os.rename(old_filtered_map_file, filtered_map_file)
        
        # Copy intermediate files to storage for index correction
        map_index_file = af.copy_to_s3(stage, id, debug, storage, BUCKET_NAME, map_index_file, True, 'map_index_files/')
        filtered_map_file = af.copy_to_s3(stage, id, debug, storage, BUCKET_NAME, filtered_map_file, True, 'filtered_map_files/')
        map_index_file = map_index_file.replace("map_index_files/", "")
        filtered_map_file = filtered_map_file.replace("filtered_map_files/", "")
                
        return fasta_chunk, fastq_chunk, map_index_file, filtered_map_file, base_name, id


    def map_alignment2(self, id, old_id, fasta_chunk, fastq_chunk, corrected_map_index_file, filtered_map_file, base_name, storage):
        """
        Second map  function, executed after the previous map function (map_alignment1) and the index correction.
        """
        
        #Global variables imported from main file
        debug = self.debug
        BUCKET_NAME = self.BUCKET_NAME
        fasta_chunks_prefix = self.fasta_chunks_prefix
        FASTA_BUCKET = self.FASTA_BUCKET
        stage = self.stage
        tolerance = self.tolerance
        file_format = self.file_format
        
        ###################################################################
        #### RECOVER DATA FROM PREVIOUS MAP
        ###################################################################
        fasta_n=re.sub(r'^\S*split_0*(\d*)\S*fasta', r'\1', fasta_chunk)
        fastq_n=re.sub('^[\s|\S]*number\':\s(\d*),[\s|\S]*$', r"\1", str(fastq_chunk))
        af.copy_to_runtime(storage, BUCKET_NAME, 'correctedIndex/', corrected_map_index_file, stage, id, debug)
        af.copy_to_runtime(storage, BUCKET_NAME, 'filtered_map_files/', filtered_map_file, stage, id, debug)
        corrected_map_index_file = "/tmp/" + corrected_map_index_file
        filtered_map_file = "/tmp/" + filtered_map_file
        fasta_chunk_folder_file = fasta_chunk.split("/")
        fasta = af.copy_to_runtime(storage, FASTA_BUCKET, fasta_chunk_folder_file[0]+"/", fasta_chunk_folder_file[1], stage, id, debug) # Download fasta chunk
        
        ###################################################################
        #### FILTER ALIGNMENTS (CORRECT .map FILE)
        ###################################################################
        map_filtering_output = sp.run(['/function/bin/map_file_index_correction.sh', corrected_map_index_file, filtered_map_file, str(tolerance)], capture_output=True)  # change to _v3.sh and runtime 20
        corrected_map_file = base_name + "_" + str(old_id) + "_filt_wline_no_corrected.map"

        ###################################################################
        #### GENERATE MPILEUP FROM MAP FILE
        ###################################################################
        mpileup_out = sp.run(['/function/bin/gempileup_run.sh', corrected_map_file, fasta], capture_output=True)
        mpileup_file = corrected_map_file + ".mpileup"

        ###################################################################
        #### FIX MPILEUP COORDINATES (pending on new fasta file splitting implementation)
        ###################################################################

        ###################################################################
        #### CONVERT MPILEUP TO PARQUET / CSV
        ###################################################################
        format_key, text_key = self.mpileup_conversion(mpileup_file, fasta_chunk, fasta_chunks_prefix, fastq_chunk, file_format, BUCKET_NAME, stage, debug, storage)
        
        return format_key, text_key
    
    
    def mpileup_conversion(self, mpileup_file, fasta_chunk, fasta_key, fastq_chunk, file_format, BUCKET_NAME, stage, debug, storage):
        # Filter mpileup file
        with open(mpileup_file, 'r') as f:
            rows = f.read().splitlines()
            content = [row.split("\t") for row in rows]
            content.pop(-1) 
            del rows

        # Convert mpileup to Pandas dataframe
        df = pd.DataFrame(data=content)
        df.columns = df.columns.astype(str)
        df['1'] = df['1'].astype(int64)

        # Remove disallowed characters
        disallowed_characters = "._-!/·¡"
        for character in disallowed_characters:
            fasta_key = fasta_key.replace(character,"")
        
        # Create intermediate key
        fasta_chunk = fasta_chunk.split("_")
        max_index = df.iloc[-1]['1']
        intermediate_key = file_format + "/" + fasta_key + "_" + fasta_chunk[-1] + "-" + fastq_chunk[0] + "_chunk" + str(fastq_chunk[1]["number"]) + "_" + str(max_index)
        
        range_index = []
        x = 0
        while x < max_index:
            if(x < max_index):
                range_index.append(x)
            x = x + 100000

        x = x + (max_index - x)
        range_index.append(x)

        df3 = df.groupby(pd.cut(df['1'], range_index)).count()
        content = ""
        for i in range(len(df3)):
            content = content + str(range_index[i+1]) + ":" + str(df3.iloc[i,0]) + "\n"

        # Write the mpileup file to the tmp directory
        map_output_file =""
        if file_format=="csv":
            map_output_file="csv/"+mpileup_file+".csv"
            df.to_csv(mpileup_file+".csv", index=False, header=False)
            #Upload the file to the s3
            with open(mpileup_file+".csv", 'rb') as f:
                storage.put_object(bucket=BUCKET_NAME, key=intermediate_key + ".csv", body=f)

            storage.put_object(bucket=BUCKET_NAME, key=intermediate_key + ".txt", body=content)
        elif file_format=="parquet":
            map_output_file="parquet/"+mpileup_file+".parquet"
            df.to_parquet(mpileup_file+".parquet")
            #Upload the file to the s3
            with open(mpileup_file+".parquet", 'rb') as f:
                storage.put_object(bucket=BUCKET_NAME, key=intermediate_key + ".parquet", body=f)
            storage.put_object(bucket=BUCKET_NAME, key=intermediate_key + ".txt", body=content)
        else:
            af.printl("file format not supported: "+file_format,stage,id,debug)
        
        return [intermediate_key+"."+file_format, intermediate_key+".txt"]
        


    ###################################################################
    ###################################################################
    #### REDUCE FUNCTIONS
    ###################################################################
    ###################################################################

    def create_sinple_key():
        import varcall_lithops_demo_v8 as config
        
        #Global variables imported from main file
        seq_type = config.seq_type
        iterdata_n = config.iterdata_n
        fastq_list = config.fastq_list
        fasta_list = config.fasta_list
        fasta_chunk_size = config.fasta_chunk_size
        fastq_chunk_size = config.fastq_chunk_size
        function_n = config.function_n
        
        
        fasta_chunk_n = math.ceil(int(iterdata_n) / len(fastq_list))
        keys = []
        fastq_file = fastq_list[0][0]

        i = 0
        for fasta_split in fasta_list:
            if i == fasta_chunk_n:
                break

            fasta_split = fasta_split.split("/")
            fasta_split = fasta_split[1]
            disallowed_characters = "._-!/·¡"
            for character in disallowed_characters:
                fasta_split = fasta_split.replace(character,"")

            fasta_split = fasta_split.split("split")
            fasta_file = fasta_split[0]
            fasta_split = fasta_split[1].split("fasta")

            fasta_file = fasta_file + "split_" + fasta_split[0] + ".fasta"
            keys.append("multipart/" + fastq_file + "-" + fasta_file + "-" + seq_type + "_" +str(fasta_chunk_size) + "fa_" + str(fastq_chunk_size) + "fq_" + function_n + ".sinple")
            i += 1

        fasta_file = fasta_file.split("split")

        keys.append("multipart/" + fastq_file + "-" + fasta_file[0]+ ".fasta" + "-" + seq_type + "_" +str(fasta_chunk_size) + "fa_" + str(fastq_chunk_size) + "fq_" + function_n + ".sinple")
        return keys


    def reduce_function(key, range, mpu_id, n_part, mpu_key, file_format, buffer_size, storage, id):
        import varcall_lithops_demo_v8 as config
        
        #Global variables imported from main file
        debug = config.debug
        BUCKET_NAME = config.BUCKET_NAME
        FASTA_BUCKET = config.FASTA_BUCKET
        
        
        stage="C"
        # check / prepare /tmp folder
        #af.file_and_folder_size("/tmp","at map function launch", stage, id, debug)
        # remove files from /tmp folder if they were left from previous runs of same script
        af.clear_tmp(stage, id, debug)
        if(os.path.exists('/tmp/reduce.mpileup')):
            os.remove('/tmp/reduce.mpileup')

        temp_mpileup = '/tmp/reduce.mpileup'
        fasta_n_first=re.sub('^\S*split_[0]*(\d*)\S*$', r"\1", key[0])
        af.printl("R\t"+fasta_n_first+"\t"+str(range['start']),stage,id,debug)

        tr0 = af.execution_time("retrieve mpileups","time_start",stage,id,debug)
        # Sequential
        # ---------------------------
        s3 = boto3.client('s3', endpoint_url="http://192.168.2.3:9000", aws_access_key_id="CR8IPBI0Q7JF7MKTF60D", aws_secret_access_key="NC21M3+4izSXLKR1WWTvl3JGJ8xT6Gwn8m0GWKqn")

        if file_format == "csv":
            expression = "SELECT * FROM s3object s WHERE cast(s._2 as int) BETWEEN %s AND %s" % (range['start'], range['end'])
            input_serialization = {'CSV': {}, 'CompressionType': 'NONE'}

        elif file_format == "parquet":
            expression = "SELECT * FROM s3object s WHERE s.\"1\" BETWEEN %s AND %s" % (range['start'], range['end'])
            input_serialization = {'Parquet': {}, 'CompressionType': 'NONE'}

        else:
            return "ERROR: Invalid format"
        
        for k in key:
            af.printl("s3 select key: "+str(k),stage,id,debug)
            fasta_n=re.sub('^\S*split_[0]*(\d*)\S*$', r"\1", k)
            fastq_n=re.sub('^\S*chunk(\d*)\S*$', r"\1", k)
            af.printl("R\t"+str(fasta_n)+"\t"+str(range['start'])+"\t"+fastq_n, stage,id,debug)
            #af.printl("fasta x fastq "+str(fasta_n)+"\t"+fastq_n, stage, id, debug)
            resp = s3.select_object_content(
                Bucket=BUCKET_NAME,
                Key=k,
                ExpressionType='SQL',
                Expression=expression,
                InputSerialization = input_serialization,
                OutputSerialization = {'CSV': {"FieldDelimiter" : "\t"}}
            )

            data = ""
            #record_count=0
            for event in resp['Payload']:
                if 'Records' in event:
                    #record_count+=1
                    records = event['Records']['Payload'].decode("UTF-8")
                    #if record_count < 10:
                        #af.printl("record "+str(record_count) + " length: "+str(len(records)), stage, id, debug)
                    data = data + records
            #af.printl("total number of mpileup records: "+str(record_count))

            #af.printl(data)

            wd = os.getcwd()
            os.chdir("/tmp")
            with open('/tmp/reduce.mpileup', 'a') as f:
                f.write(data)
            os.chdir(wd)
            del data
        # ---------------------------
        tr1 = af.execution_time("retrieve mpileups","time_end",stage,id,debug)
        af.printl(f' retrieve mpileups: execution_time:  {tr1 - tr0}: s',stage,id,debug)
        af.printl("r\t"+fasta_n_first+"\t"+str(range['start']),stage,id,debug)

        chr_table = af.copy_to_runtime(storage, FASTA_BUCKET, "cache/", "hashtable.data", stage, id, debug)
        af.printl("chromosome hash table: "+str(chr_table), stage, id, debug)
        with open(chr_table, 'r') as f:
            for line in f:
                af.printl(line, stage, id, debug)

        wd = os.getcwd()
        os.chdir("/tmp")
        
        tr2 = af.execution_time("merge mpileups","time_start",stage,id,debug)
        af.printl("S\t"+fasta_n_first+"\t"+str(range['start']),stage,id,debug)
        #af.file_and_folder_size("/tmp","before mpileup_merge_reducev3.sh", stage, id, debug)
        #af.copy_to_s3(stage, id, debug, storage, BUCKET_NAME, temp_mpileup, temp_to_s3, s3_temp_folder)
        af.printl("Starting Merge Script",stage,id,debug)
        sinple_out = sp.check_output(['bash', '/function/bin/mpileup_merge_reducev3_nosinple.sh', temp_mpileup, '/function/bin', buffer_size])
        sinple_out = sinple_out.decode('UTF-8')
        af.printl("Finished Merge Script",stage,id,debug)
        sinple_name=temp_mpileup+'_merged.mpileup'

        # write output to /tmp
        with open(sinple_name, 'w') as f:
            f.write(sinple_out)
        tr3 = af.execution_time("merge mpileups","time_end",stage,id,debug)
        af.printl("s\t"+fasta_n_first+"\t"+str(range['start']),stage,id,debug)
        af.printl(f' merge mpileups: execution_time:  {tr3 - tr2}: s',stage,id,debug)
        
        #af.file_and_folder_size("/tmp","after mpileup_merge_reducev3.sh", stage, id, debug)
        # copy output to s3
        #af.copy_to_s3(stage, id, debug, storage, BUCKET_NAME, sinple_name, temp_to_s3, s3_temp_folder)

        os.chdir(wd)
        #Upload part
        s3 = boto3.client('s3', endpoint_url="http://192.168.2.3:9000", aws_access_key_id="CR8IPBI0Q7JF7MKTF60D", aws_secret_access_key="NC21M3+4izSXLKR1WWTvl3JGJ8xT6Gwn8m0GWKqn")
        part = s3.upload_part(
            Body = sinple_out,
            Bucket = BUCKET_NAME,
            Key = mpu_key,
            UploadId = mpu_id,
            PartNumber = n_part
        )
        tr4 = time.time()
        af.printl(f'reduce 2 - lambda: execution_time_total: {tr4 - tr0}: s',stage,id,debug)
        return {"PartNumber" : n_part, "ETag" : part["ETag"], "mpu_id": mpu_id}