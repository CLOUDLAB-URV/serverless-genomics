import os
import re
import multiprocessing
import subprocess as sp
from typing import Tuple
import pandas as pd
from numpy import int64
import fastq_functions as fq_func
import aux_functions as aux
from varcall_arguments import Arguments
from lithops import Storage
import zipfile


class AlignmentMapper:
    def __init__(self, fasta_chunks_prefix: str, args: Arguments):
        self.fasta_chunks_prefix = fasta_chunks_prefix
        self.args = args
    
    ###################################################################
    ###################################################################
    #### MAP FUNCTIONS
    ###################################################################
    ###################################################################


    def map_alignment1(self, id: int, fasta_chunk: dict, fastq_chunk: str, storage: Storage):
        """
        First map function to filter and map fasta + fastq chunks. Some intermediate files
        are uploaded into the cloud storage for subsequent index correction, after which
        the final part of the map function (map_alignment2) can be executed.
        """
        
        # CONTROL VARIABLES 
        fastq_n=re.sub('^[\s|\S]*number\':\s(\d*),[\s|\S]*$', r"\1", str(fastq_chunk))
        func_key=self.args.fq_seqname+"_fq"+fastq_n+"-fa"+str(fasta_chunk['id'])

        
        ###################################################################
        #### PROCESSING FASTQ CHUNKS
        ###################################################################
        # Download fastq chunk depending on source
        fastq1, fastq2, base_name = self.download_fastq(fastq_chunk)
        
        ###################################################################
        #### PROCESSING FASTA CHUNKS
        ###################################################################
        fasta_folder_file = fasta_chunk['key_fasta'].split("/")
        fasta = aux.copy_to_runtime(storage, self.args.fasta_bucket, fasta_folder_file[0]+"/", fasta_folder_file[1], 
                                {'Range': f"bytes={fasta_chunk['chunk'][0]['offset_base']}-{fasta_chunk['chunk'][1]['last_byte+']}"}, fasta_chunk) # Download fasta chunk 
    
        gem_ref_nosuffix = os.path.splitext(fasta)[0]
        gem_ref = gem_ref_nosuffix + '.gem'
        cpus=multiprocessing.cpu_count()
        sp.run(['gem-indexer', '--input', fasta, '--threads', str(cpus), '-o', gem_ref_nosuffix], capture_output=True)

        ###################################################################
        #### GENERATE ALIGNMENT AND ALIGNMENT INDEX (FASTQ TO MAP)
        ###################################################################
        sp.run(['/function/bin/map_index_and_filter_map_file_cmd_awsruntime.sh', gem_ref, fastq1, fastq2, base_name, self.args.datasource, self.args.seq_type], capture_output=True)

        # Reorganize file names
        map_index_file = base_name + "_map.index.txt"
        os.rename(map_index_file,"/tmp/" + func_key+ "_map.index.txt")
        map_index_file = "/tmp/" + func_key+ "_map.index.txt"
        old_filtered_map_file = base_name + "_filt_wline_no.map"
        filtered_map_file = base_name + "_" + str(id) + "_filt_wline_no.map"
        os.rename(old_filtered_map_file, filtered_map_file)
        

        # Compress the filtered map file
        zipname = filtered_map_file + ".zip"
        with zipfile.ZipFile(zipname, 'w', compression=zipfile.ZIP_BZIP2, compresslevel=9) as zf:
            zf.write(filtered_map_file)
        
        # Copy intermediate files to storage for index correction       
        filtered_map_file = aux.copy_to_s3(storage, self.args.bucket, zipname, True, 'filtered_map_files/')
        filtered_map_file = filtered_map_file.replace("filtered_map_files/", "")
        
        map_index_file = aux.copy_to_s3(storage, self.args.bucket, map_index_file, True, 'map_index_files/')
        map_index_file = map_index_file.replace("map_index_files/", "")
                
        return fasta_chunk, fastq_chunk, map_index_file, filtered_map_file, base_name, id


    def map_alignment2(self, old_id: int, fasta_chunk: dict, fastq_chunk: str, corrected_map_index_file: str, filtered_map_file: str, base_name: str, storage: Storage):

        """
        Second map  function, executed after the previous map function (map_alignment1) and the index correction.
        """
        
        ###################################################################
        #### RECOVER DATA FROM PREVIOUS MAP
        ###################################################################
        corrected_map_index_file = aux.copy_to_runtime(storage, self.args.bucket, 'corrected_index/', corrected_map_index_file)
        filtered_map_file = aux.copy_to_runtime(storage, self.args.bucket, 'filtered_map_files/', filtered_map_file)
        
        fasta_folder_file = fasta_chunk['key_fasta'].split("/") 
        fasta = aux.copy_to_runtime(storage, self.args.fasta_bucket, fasta_folder_file[0]+"/", fasta_folder_file[1], 
                                {'Range': f"bytes={fasta_chunk['chunk'][0]['offset_base']}-{fasta_chunk['chunk'][1]['last_byte+']}"}, fasta_chunk) # Download fasta chunk
        
        with zipfile.ZipFile(filtered_map_file) as zf:
            zf.extractall("/")
        filtered_map_file = filtered_map_file.replace(".zip", "")
        
        ###################################################################
        #### FILTER ALIGNMENTS (CORRECT .map FILE)
        ###################################################################
        map_filtering_output = sp.run(['/function/bin/map_file_index_correction.sh', corrected_map_index_file, filtered_map_file, str(self.args.tolerance)], capture_output=True)  # change to _v3.sh and runtime 20
        corrected_map_file = base_name + "_" + str(old_id) + "_filt_wline_no_corrected.map"

        ###################################################################
        #### GENERATE MPILEUP FROM MAP FILE
        ###################################################################
        mpileup_out = sp.run(['/function/bin/gempileup_run.sh', corrected_map_file, fasta], capture_output=True)
        mpileup_file = corrected_map_file + ".mpileup"

        ###################################################################
        #### CONVERT MPILEUP TO PARQUET / CSV
        ###################################################################
        format_key, text_key = self.mpileup_conversion(mpileup_file, fasta_chunk, fastq_chunk, storage)
        
        return format_key, text_key
        

    ###################################################################
    ###################################################################
    #### AUXILIARY MAP FUNCTIONS
    ###################################################################
    ###################################################################
    
    def download_fastq(self, fastq_chunk) -> Tuple[str, str, str]:
        """
        Download fastq chunks depending on source
        """   
        if self.args.seq_type == "paired-end":
            fastq1 = fq_func.fastq_to_mapfun(fastq_chunk[0], fastq_chunk[1])
            base_name = os.path.splitext(fastq1)[0] #+'.pe'
            fastq2 = "yes"
        else:   # single-end sequencing
            fastq1 = fq_func.fastq_to_mapfun(fastq_chunk[0], fastq_chunk[1])
            fastq2 = "no"
            base_name = os.path.splitext(fastq1)[0] #+'.se'
        return fastq1, fastq2, base_name
    
    
    def mpileup_conversion(self, mpileup_file: str, fasta_chunk: dict, fastq_chunk: str, storage: Storage) -> Tuple [str, str]:

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
        fasta_key = self.fasta_chunks_prefix
        disallowed_characters = "._-!/·¡"
        for character in disallowed_characters:
            fasta_key = fasta_key.replace(character,"")
        
        # Create intermediate key
        fasta_chunk = str(fasta_chunk['id'])
        max_index = df.iloc[-1]['1']
        intermediate_key = self.args.file_format + "/" + fasta_key + "_" + fasta_chunk + "-" + fastq_chunk[0] + "_chunk" + str(fastq_chunk[1]["number"]) + "_" + str(max_index)

        
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
        storage.put_object(bucket=self.args.bucket, key=intermediate_key + ".txt", body=content)
        
        # Write the mpileup file to the tmp directory
        if self.args.file_format=="csv":
            df.to_csv(mpileup_file+".csv", index=False, header=False)
            with open(mpileup_file+".csv", 'rb') as f:
                storage.put_object(bucket=self.args.bucket, key=intermediate_key + ".csv", body=f)
            
        elif self.args.file_format=="parquet":
            df.to_parquet(mpileup_file+".parquet")
            with open(mpileup_file+".parquet", 'rb') as f:
                storage.put_object(bucket=self.args.bucket, key=intermediate_key + ".parquet", body=f)
        
        return [intermediate_key+"."+self.args.file_format, intermediate_key+".txt"]