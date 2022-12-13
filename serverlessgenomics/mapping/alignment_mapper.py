import bz2
import logging
import os
import multiprocessing
import pathlib
import shutil
import subprocess as sp
import tempfile
from typing import Tuple
import pandas as pd
from numpy import int64
from pathlib import PurePosixPath

from .data_fetch import fetch_fastq_chunk, fetch_fasta_chunk
from ..utils import copy_to_runtime, force_delete_local_path
from ..parameters import PipelineRun
from lithops import Storage
from lithops.storage.utils import StorageNoSuchKeyError
import zipfile

logger = logging.getLogger(__name__)


def gem_indexer_mapper(pipeline_params: PipelineRun,
                       fasta_chunk_id: int, fasta_chunk: dict,
                       fastq_chunk_id: int, fastq_chunk: dict,
                       storage: Storage):
    """
    Lithops callee function
    First map function to filter and map fasta + fastq chunks. Some intermediate files
    are uploaded into the cloud storage for subsequent index correction, after which
    the final part of the map function (map_alignment2) can be executed.

    Args:
        mapper_id (int): mapper id
        fasta_chunk (dict): contains all the necessary info related to the fasta chunk
        fastq_chunk (str): fastq chunk key
        exec_param (str): string used to differentiate this pipeline execution from others with different parameters
        storage (Storage): s3 storage instance, generated by lithops

    Returns:
        Tuple[str]: multiple values needed for index correction and the second map phase
    """

    mapper_id = f'fa{fasta_chunk_id}-fq{fastq_chunk_id}'
    base_name = 'SRRXXXXXX'
    map_index_key = os.path.join(pipeline_params.tmp_prefix, pipeline_params.run_id, 'gem-mapper',
                                 f'fq{fastq_chunk_id}', f'fa{fasta_chunk_id}', base_name + '_map.index.txt.bz2')
    filtered_map_key = os.path.join(pipeline_params.tmp_prefix, pipeline_params.run_id, 'gem-mapper',
                                    f'fq{fastq_chunk_id}', f'fa{fasta_chunk_id}',
                                    base_name + '_filt_wline_no.map.bz2')

    # Check if output files already exist in storage
    try:
        storage.head_object(bucket=pipeline_params.storage_bucket, key=map_index_key)
        storage.head_object(bucket=pipeline_params.storage_bucket, key=filtered_map_key)
        # If they exist, return the keys and skip computing this chunk
        return fastq_chunk_id, fasta_chunk_id, map_index_key, filtered_map_key
    except StorageNoSuchKeyError:
        # If any output is missing, proceed
        pass

    # Make temp dir and ch into it, save pwd to restore it later
    tmp_dir = tempfile.mkdtemp()
    pwd = os.getcwd()
    os.chdir(tmp_dir)

    try:
        # Get fastq chunk and store it to disk in tmp directory
        fastq_chunk_filename = f"chunk_{fastq_chunk['chunk_id']}.fastq"
        fastqgz_idx_key, _ = pipeline_params.fastqgz_idx_keys
        fetch_fastq_chunk(fastq_chunk, fastq_chunk_filename, storage, pipeline_params.fastq_path,
                          pipeline_params.storage_bucket, fastqgz_idx_key)

        # Get fasta chunk and store it to disk in tmp directory
        fasta_chunk_filename = f"chunk_{fasta_chunk['chunk_id']}.fasta"
        fetch_fasta_chunk(fasta_chunk, fasta_chunk_filename, storage, pipeline_params.fasta_path)

        # TODO try to move fasta indexing to another preprocessing step
        gem_index_filename = os.path.join(f'{mapper_id}.gem')
        cpus = multiprocessing.cpu_count()

        # gem-indexer appends .gem to output file
        cmd = ['gem-indexer', '--input', fasta_chunk_filename, '--threads', str(cpus), '-o',
               gem_index_filename.replace('.gem', '')]
        print(' '.join(cmd))
        out = sp.run(cmd, capture_output=True)
        print(out.stderr.decode('utf-8'))
        # TODO apparently 1 is return code for success (why)
        assert out.returncode == 1

        # GENERATE ALIGNMENT AND ALIGNMENT INDEX (FASTQ TO MAP)
        # TODO refactor bash script
        # TODO support implement paired-end, replace not-used with 2nd fastq chunk
        # TODO use proper tmp directory instead of uuid base name
        # TODO add support for sra source
        cmd = ['/function/bin/map_index_and_filter_map_file_cmd_awsruntime.sh', gem_index_filename,
               fastq_chunk_filename, "not-used", base_name, "s3", "single-end"]
        print(' '.join(cmd))
        out = sp.run(cmd, capture_output=True)
        print(out.stdout.decode('utf-8'))
        print(out.stderr.decode('utf-8'))

        # Reorganize file names
        map_index_filename = os.path.join(tmp_dir, base_name + "_map.index.txt")
        shutil.move(base_name + "_map.index.txt", map_index_filename)

        filtered_map_filename = os.path.join(tmp_dir, base_name + "_" + str(mapper_id) + "_filt_wline_no.map")
        shutil.move(base_name + "_filt_wline_no.map", filtered_map_filename)

        # Compress outputs
        zipped_map_index_filename = map_index_filename + ".bz2"
        with zipfile.ZipFile(zipped_map_index_filename, 'w', compression=zipfile.ZIP_BZIP2, compresslevel=9) as zf:
            zf.write(map_index_filename, arcname=PurePosixPath(map_index_filename).name)

        zipped_filtered_map_filename = filtered_map_filename + ".bz2"
        with zipfile.ZipFile(zipped_filtered_map_filename, 'w', compression=zipfile.ZIP_BZIP2, compresslevel=9) as zf:
            zf.write(filtered_map_filename, arcname=PurePosixPath(filtered_map_filename).name)

        # Copy result files to storage for index correction
        storage.upload_file(file_name=zipped_map_index_filename, bucket=pipeline_params.storage_bucket,
                            key=map_index_key)
        storage.upload_file(file_name=zipped_filtered_map_filename, bucket=pipeline_params.storage_bucket,
                            key=filtered_map_key)
    finally:
        os.chdir(pwd)
        force_delete_local_path(tmp_dir)

    return fastq_chunk_id, fasta_chunk_id, map_index_key, filtered_map_key


def index_correction(pipeline_params, fastq_chunk_id, map_index_keys, storage: Storage):
    """
    Lithops callee function
    Corrects the index after the first map iteration.
    All the set files must have the prefix "map_index_files/".
    Corrected indices will be stored with the prefix "corrected_index/".

    Args:
        setname (str): files to be corrected
        bucket (str): s3 bucket where the set is stored
        exec_param (str): string used to differentiate this pipeline execution from others with different parameters
        storage (Storage): s3 storage instance, generated by lithops
    """
    set_name = f'fq_{fastq_chunk_id}'
    # TODO get base name from params
    base_name = 'SRRXXXXXX'
    pwd = os.getcwd()

    # TODO replace with a proper file name (maybe indicating fastq chunk id)
    output_file = 'merged_filtered_index.txt'
    zipped_output_file = output_file + ".bz2"
    corrected_index_key = os.path.join(pipeline_params.tmp_prefix, pipeline_params.run_id, 'index-correction',
                                       f'fq{fastq_chunk_id}', zipped_output_file)

    # Check if output files already exist in storage
    try:
        storage.head_object(bucket=pipeline_params.storage_bucket, key=corrected_index_key)
        # If they exist, return the keys and skip computing this chunk
        return fastq_chunk_id, corrected_index_key
    except StorageNoSuchKeyError:
        # If the output is missing, proceed
        pass

    # Download all map files for this fastq chunk
    input_temp_dir = tempfile.mkdtemp()
    output_temp_dir = tempfile.mkdtemp()
    try:
        for i, map_index_key in enumerate(map_index_keys):
            local_compressed_map_path = os.path.join(input_temp_dir, f'map_{i}.map.bz2')
            storage.download_file(bucket=pipeline_params.storage_bucket, key=map_index_key, file_name=local_compressed_map_path)
            with zipfile.ZipFile(local_compressed_map_path, 'r', compression=zipfile.ZIP_BZIP2, compresslevel=9) as zf:
                zf.extractall(input_temp_dir)
            # TODO set proper map index file name
            os.rename(os.path.join(input_temp_dir, f'{base_name}_map.index.txt'), os.path.join(input_temp_dir, f'{i}_map.index.txt'))
            os.remove(local_compressed_map_path)
        print(os.listdir(input_temp_dir))

        # TODO chdir required or binary_reducer.sh
        os.chdir(output_temp_dir)

        # Execute correction scripts
        intermediate_file = f'{set_name}.intermediate.txt'
        # TODO binary_reducer.sh could be more efficiently implemented using python and consuming files from storage as a generator
        cmd = f'/function/bin/binary_reducer.sh /function/bin/merge_gem_alignment_metrics.sh 4 {input_temp_dir}/* > {intermediate_file}'
        print(cmd)
        proc = sp.run(cmd, shell=True, universal_newlines=True, capture_output=True)
        print(proc.stdout)
        print(proc.stderr)
        proc.check_returncode()

        cmd = f'/function/bin/filter_merged_index.sh {intermediate_file} {output_file}'
        print(cmd)
        proc = sp.run(cmd, shell=True, check=True, universal_newlines=True)
        print(proc.stdout)
        print(proc.stderr)
        proc.check_returncode()

        print(os.listdir())

        # Compress output
        with zipfile.ZipFile(zipped_output_file, 'w', compression=zipfile.ZIP_BZIP2, compresslevel=9) as zf:
            zf.write(output_file, arcname=output_file)

        # Upload corrected index to storage
        storage.upload_file(bucket=pipeline_params.storage_bucket, key=corrected_index_key, file_name=zipped_output_file)

        os.chdir(pwd)
        return fastq_chunk_id, corrected_index_key
    finally:
        os.chdir(pwd)
        force_delete_local_path(input_temp_dir)
        force_delete_local_path(output_temp_dir)


def filter_index_to_mpileup(pipeline_params, fasta_chunk_id, fasta_chunk, fastq_chunk_id, fastq_chunk,
                            filtered_map_key, corrected_index_key, storage):
    """
    Second map  function, executed after the previous map function (map_alignment1) and the index correction.

    Args:
        old_id (int): id used in the previous map function
        fasta_chunk (dict): contains all the necessary info related to the fasta chunk
        fastq_chunk (str): fastq chunk key
        corrected_map_index_file (str): key of the corrected index stored in s3
        filtered_map_file (str): key of the filtered map file generated in the previous map function and stored in s3
        base_name (str): name used in the name of some of the generated files
        exec_param (str): string used to differentiate this pipeline execution from others with different parameters
        storage (Storage): s3 storage instance, generated by lithops

    Returns:
        Tuple[str]: keys to the generated txt and csv/parquet files (stored in s3)
    """
    temp_dir = tempfile.mkdtemp()
    pwd = os.getcwd()
    # TODO get base name from params
    base_name = 'SRRXXXXXX'

    corrected_map_file = f'{base_name}_fa{fasta_chunk_id}-fq{fastq_chunk_id}_filt_wline_no_corrected.map'
    mpileup_file = corrected_map_file + ".mpileup"
    mpipleup_key = os.path.join(pipeline_params.tmp_prefix, pipeline_params.run_id, 'mpileups',
                                f'fq{fastq_chunk_id}', f'fa{fasta_chunk_id}', mpileup_file)
    
    # Check if output files already exist in storage
    try:
        storage.head_object(bucket=pipeline_params.storage_bucket, key=mpipleup_key)
        # If they exist, return the keys and skip computing this chunk
        return fastq_chunk_id, mpipleup_key
    except StorageNoSuchKeyError:
        # If the output is missing, proceed
        pass

    try:
        os.chdir(temp_dir)

        # Recover fasta chunk
        fasta_chunk_filename = f"chunk_{fasta_chunk['chunk_id']}.fasta"
        fetch_fasta_chunk(fasta_chunk, fasta_chunk_filename, storage, pipeline_params.fasta_path)

        # Recover filtered map file
        bz2_filt_map_filename = pathlib.PurePosixPath(filtered_map_key).name
        storage.download_file(bucket=pipeline_params.storage_bucket, key=filtered_map_key,
                              file_name=bz2_filt_map_filename)
        with zipfile.ZipFile(bz2_filt_map_filename, 'r', compression=zipfile.ZIP_BZIP2, compresslevel=9) as zf:
            zf.extractall()
        os.remove(bz2_filt_map_filename)

        # Get corrected map index for this fastq chunk
        bz2_corrected_index_filename = pathlib.PurePosixPath(corrected_index_key).name
        storage.download_file(bucket=pipeline_params.storage_bucket, key=corrected_index_key,
                              file_name=bz2_corrected_index_filename)
        with zipfile.ZipFile(bz2_corrected_index_filename, 'r', compression=zipfile.ZIP_BZIP2, compresslevel=9) as zf:
            zf.extractall()
        os.remove(bz2_corrected_index_filename)

        print(os.listdir(temp_dir))

        # Filter aligments with corrected map file
        filt_map_filename = f'{base_name}_fa{fasta_chunk_id}-fq{fastq_chunk_id}_filt_wline_no.map'
        # TODO replace with a proper file name (maybe indicating fastq chunk id)
        corrected_index_filename = 'merged_filtered_index.txt'

        cmd = ['/function/bin/map_file_index_correction.sh', corrected_index_filename,
               filt_map_filename, str(pipeline_params.tolerance)]
        print(' '.join(cmd))
        proc = sp.run(cmd, capture_output=True)  # change to _v3.sh and runtime 20
        print(proc.stdout.decode('utf-8'))
        print(proc.stderr.decode('utf-8'))

        print(os.listdir(temp_dir))

        # Generate mpileup
        cmd = ['/function/bin/gempileup_run.sh', corrected_map_file, fasta_chunk_filename]
        print(' '.join(cmd))
        proc = sp.run(cmd, capture_output=True)
        print(proc.stdout.decode('utf-8'))
        print(proc.stderr.decode('utf-8'))

        # Store output to storage
        mpipleup_key = os.path.join(pipeline_params.tmp_prefix, pipeline_params.run_id, 'mpileups',
                                    f'fq{fastq_chunk_id}', f'fa{fasta_chunk_id}', mpileup_file)
        storage.upload_file(bucket=pipeline_params.storage_bucket, key=corrected_index_key, file_name=mpipleup_key)

        return fastq_chunk_id, fasta_chunk_id, mpipleup_key
    finally:
        os.chdir(pwd)
        force_delete_local_path(temp_dir)


def mpileup_conversion(self, mpileup_file: str, fasta_chunk: dict, fastq_chunk: str, exec_param: str,
                       storage: Storage) -> Tuple[str]:
    """
    Convert resulting data to csv/parquet and txt

    Args:
        mpileup_file (str): mpileup file
        fasta_chunk (dict): contains all the necessary info related to the fasta chunk
        fastq_chunk (str): fastq chunk key
        exec_param (str): string used to differentiate this pipeline execution from others with different parameters
        storage (Storage): s3 storage instance, generated by lithops

    Returns:
        Tuple[str]: keys to the generated txt and csv/parquet files (stored in s3)
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
        fasta_key = fasta_key.replace(character, "")

    # Create intermediate key
    fasta_chunk = str(fasta_chunk['id'])
    max_index = df.iloc[-1]['1']
    intermediate_key = self.args.file_format + "/" + exec_param + "/" + fasta_key + "_" + fasta_chunk + "-" + \
                       fastq_chunk[0] + "_chunk" + str(fastq_chunk[1]["number"]) + "_" + str(max_index)

    range_index = []
    x = 0
    while x < max_index:
        if (x < max_index):
            range_index.append(x)
        x = x + 100000

    x = x + (max_index - x)
    range_index.append(x)

    df3 = df.groupby(pd.cut(df['1'], range_index)).count()
    content = ""
    for i in range(len(df3)):
        content = content + str(range_index[i + 1]) + ":" + str(df3.iloc[i, 0]) + "\n"

    # Upload .txt file to storage
    storage.put_object(bucket=self.args.storage_bucket, key=intermediate_key + ".txt", body=content)

    # Write the mpileup file to the tmp directory
    if self.args.file_format == "csv":
        df.to_csv(mpileup_file + ".csv", index=False, header=False)
        with open(mpileup_file + ".csv", 'rb') as f:
            storage.put_object(bucket=self.args.storage_bucket, key=intermediate_key + ".csv", body=f)

    elif self.args.file_format == "parquet":
        df.to_parquet(mpileup_file + ".parquet")
        with open(mpileup_file + ".parquet", 'rb') as f:
            storage.put_object(bucket=self.args.storage_bucket, key=intermediate_key + ".parquet", body=f)

    return [intermediate_key + "." + self.args.file_format, intermediate_key + ".txt"]
