import os.path
import re
from random import randint
import multiprocessing
import subprocess as sp
import pandas as pd
from numpy import int64
import lithopsgenetics
import lithopsgenetics.auxiliaryfunctions as af


class MapFunctions:
    def __init__(self, fq_seqname, datasource, seq_type, debug, BUCKET_NAME, fastq_folder, idx_folder,
                  fasta_chunks_prefix, FASTA_BUCKET, stage, file_format, tolerance):
        self.fq_seqname = fq_seqname
        self.datasource = datasource
        self.seq_type = seq_type
        self.debug = debug
        self.BUCKET_NAME = BUCKET_NAME
        self.fastq_folder = fastq_folder
        self.idx_folder = idx_folder
        self.fasta_chunks_prefix = fasta_chunks_prefix
        self.FASTA_BUCKET = FASTA_BUCKET
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
        the final part of the map function (map_alignment2) can be executed.
        """
        
        # CONTROL VARIABLES 
        stage="A"+str(randint(1000,9999))
        fasta_n=re.sub(r'^\S*split_0*(\d*)\S*fasta', r'\1', fasta_chunk)
        fastq_n=re.sub('^[\s|\S]*number\':\s(\d*),[\s|\S]*$', r"\1", str(fastq_chunk))
        func_key=self.fq_seqname+"_fq"+fastq_n+"-fa"+fasta_n
        
        ###################################################################
        #### PROCESSING FASTQ CHUNKS
        ###################################################################
        # Download fastq chunk depending on source
        fastq1 =""
        fastq2 =""
        fastq1, fastq2, base_name = self.download_fastq(id, stage, fastq_chunk)
        
        ###################################################################
        #### PROCESSING FASTA CHUNKS
        ###################################################################
        fasta_chunk_folder_file = fasta_chunk.split("/")
        fasta = af.copy_to_runtime(storage, self.FASTA_BUCKET, fasta_chunk_folder_file[0]+"/", fasta_chunk_folder_file[1], stage, id, self.debug)
        gem_ref_nosuffix = os.path.splitext(fasta)[0]
        gem_ref = gem_ref_nosuffix + '.gem'
        cpus=multiprocessing.cpu_count()
        sp.run(['gem-indexer', '--input', fasta, '--threads', str(cpus), '-o', gem_ref_nosuffix], capture_output=True)


        ###################################################################
        #### GENERATE ALIGNMENT AND ALIGNMENT INDEX (FASTQ TO MAP)
        ###################################################################
        sp.run(['/function/bin/map_index_and_filter_map_file_cmd_awsruntime.sh', gem_ref, fastq1, fastq2, base_name, self.datasource, self.seq_type], capture_output=True)

        # Reorganize file names
        map_index_file = base_name + "_map.index.txt"
        os.rename(map_index_file,"/tmp/" + func_key+ "_map.index.txt")
        map_index_file = "/tmp/" + func_key+ "_map.index.txt"
        old_filtered_map_file = base_name + "_filt_wline_no.map"
        filtered_map_file = base_name + "_" + str(id) + "_filt_wline_no.map"
        os.rename(old_filtered_map_file, filtered_map_file)
        
        # Copy intermediate files to storage for index correction
        map_index_file = af.copy_to_s3(stage, id, self.debug, storage, self.BUCKET_NAME, map_index_file, True, 'map_index_files/')
        filtered_map_file = af.copy_to_s3(stage, id, self.debug, storage, self.BUCKET_NAME, filtered_map_file, True, 'filtered_map_files/')
        map_index_file = map_index_file.replace("map_index_files/", "")
        filtered_map_file = filtered_map_file.replace("filtered_map_files/", "")
                
        return fasta_chunk, fastq_chunk, map_index_file, filtered_map_file, base_name, id


    def map_alignment2(self, old_id, fasta_chunk, fastq_chunk, corrected_map_index_file, filtered_map_file, base_name, storage):
        """
        Second map  function, executed after the previous map function (map_alignment1) and the index correction.
        """
        
        ###################################################################
        #### RECOVER DATA FROM PREVIOUS MAP
        ###################################################################
        fasta_n=re.sub(r'^\S*split_0*(\d*)\S*fasta', r'\1', fasta_chunk)
        fastq_n=re.sub('^[\s|\S]*number\':\s(\d*),[\s|\S]*$', r"\1", str(fastq_chunk))
        af.copy_to_runtime(storage, self.BUCKET_NAME, 'correctedIndex/', corrected_map_index_file, self.stage, old_id, self.debug)
        af.copy_to_runtime(storage, self.BUCKET_NAME, 'filtered_map_files/', filtered_map_file, self.stage, old_id, self.debug)
        corrected_map_index_file = "/tmp/" + corrected_map_index_file
        filtered_map_file = "/tmp/" + filtered_map_file
        fasta_chunk_folder_file = fasta_chunk.split("/")
        fasta = af.copy_to_runtime(storage, self.FASTA_BUCKET, fasta_chunk_folder_file[0]+"/", fasta_chunk_folder_file[1], self.stage, old_id, self.debug) # Download fasta chunk
        
        ###################################################################
        #### FILTER ALIGNMENTS (CORRECT .map FILE)
        ###################################################################
        map_filtering_output = sp.run(['/function/bin/map_file_index_correction.sh', corrected_map_index_file, filtered_map_file, str(self.tolerance)], capture_output=True)  # change to _v3.sh and runtime 20
        corrected_map_file = base_name + "_" + str(old_id) + "_filt_wline_no_corrected.map"

        ###################################################################
        #### GENERATE MPILEUP FROM MAP FILE
        ###################################################################
        mpileup_out = sp.run(['/function/bin/gempileup_run.sh', corrected_map_file, fasta], capture_output=True)
        mpileup_file = corrected_map_file + ".mpileup"

        ###################################################################
        #### CONVERT MPILEUP TO PARQUET / CSV
        ###################################################################
        format_key, text_key = self.mpileup_conversion(mpileup_file, fasta_chunk, self.fasta_chunks_prefix, fastq_chunk, self.file_format, self.BUCKET_NAME, storage)
        
        return format_key, text_key
        

    ###################################################################
    ###################################################################
    #### AUXILIARY MAP FUNCTIONS
    ###################################################################
    ###################################################################
    
    def download_fastq(self, id, stage, fastq_chunk):
        """
        Download fastq chunks depending on source
        """
        if self.datasource == "s3":
            if self.seq_type == "paired-end":
                fastq1 = lithopsgenetics.fastq_to_mapfun("fastq1", fastq_chunk[0][0], fastq_chunk[0][1], self.BUCKET_NAME, self.fastq_folder, self.idx_folder, self.datasource, stage,id, self.debug)
                fastq2 = lithopsgenetics.fastq_to_mapfun("fastq2", fastq_chunk[1][0], fastq_chunk[1][1], self.BUCKET_NAME, self.fastq_folder, self.idx_folder, self.datasource, stage,id, self.debug)
                base_name = os.path.splitext(fastq1)[0]#+'.pe'
            else:   # single-end sequencing
                fastq1 = lithopsgenetics.fastq_to_mapfun("fastq", fastq_chunk[0], fastq_chunk[1], self.BUCKET_NAME, self.fastq_folder, self.idx_folder, self.datasource, stage,id, self.debug)
                base_name = os.path.splitext(fastq1)[0]+'.se'
                fastq2 = "no"      
        elif self.datasource == "SRA":
            if self.seq_type == "paired-end":
                fastq1 = lithopsgenetics.fastq_to_mapfun("fastq", fastq_chunk[0], fastq_chunk[1], self.BUCKET_NAME, self.fastq_folder, self.idx_folder, self.datasource, stage,id, self.debug)
                base_name = os.path.splitext(fastq1)[0]#+'.pe'
                fastq2 = "yes"
            else:   # single-end sequencing
                fastq1 = lithopsgenetics.fastq_to_mapfun("fastq", fastq_chunk[0], fastq_chunk[1], self.BUCKET_NAME, self.fastq_folder, self.idx_folder, self.datasource, stage,id, self.debug)
                fastq2 = "no"
                base_name = os.path.splitext(fastq1)[0]#+'.se'
        return fastq1, fastq2, base_name
    
    
    def mpileup_conversion(self, mpileup_file, fasta_chunk, fasta_key, fastq_chunk, file_format, BUCKET_NAME, storage):
        """
        Convert resulting data to csv/parquet and txt
        """
        
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
            
        # Upload .txt file to storage
        storage.put_object(bucket=BUCKET_NAME, key=intermediate_key + ".txt", body=content)
        
        # Write the mpileup file to the tmp directory
        if file_format=="csv":
            df.to_csv(mpileup_file+".csv", index=False, header=False)
            with open(mpileup_file+".csv", 'rb') as f:
                storage.put_object(bucket=BUCKET_NAME, key=intermediate_key + ".csv", body=f)
            
        elif file_format=="parquet":
            df.to_parquet(mpileup_file+".parquet")
            with open(mpileup_file+".parquet", 'rb') as f:
                storage.put_object(bucket=BUCKET_NAME, key=intermediate_key + ".parquet", body=f)
        
        return [intermediate_key+"."+file_format, intermediate_key+".txt"]