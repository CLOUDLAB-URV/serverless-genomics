import logging
import os
import multiprocessing
import pathlib
import shutil
import subprocess as sp
import tempfile
from functools import partial
from typing import Tuple
import pandas as pd
from numpy import int64
from pathlib import PurePosixPath
from time import time

from ..datasource import fetch_fasta_chunk, fetch_fastq_chunk
from ..datasource.fetch import fetch_gem_chunk
from ..utils import force_delete_local_path, get_storage_tmp_prefix
from ..pipeline import PipelineParameters
from ..stats import Stats
from lithops import Storage
from lithops.storage.utils import StorageNoSuchKeyError
import zipfile

logger = logging.getLogger(__name__)


def align_mapper(
    pipeline_params: PipelineParameters,
    run_id: str,
    mapper_id: str,
    fasta_chunk: dict,
    fastq_chunk: dict,
    storage: Storage,
):
    """
    Lithops callee function
    First map function to filter and map fasta + fastq chunks. Some intermediate files
    are uploaded into the cloud storage for subsequent index correction, after which
    the final part of the map function (map_alignment2) can be executed.
    """
    print("starting align_mapper")
    stats = Stats()
    stats.set_value("mapper_id", mapper_id)
    stats.start_timer("function")

    # tmp prefix generator for this mapper
    mapper_storage_tmp_prefix = partial(get_storage_tmp_prefix, run_id, "align_mapper", mapper_id)

    map_index_key = mapper_storage_tmp_prefix(pipeline_params.sra_accession + "_map.index.txt.bz2")
    filtered_map_key = mapper_storage_tmp_prefix(pipeline_params.sra_accession + "_filt_wline_no.map.bz2")

    # Check if output files already exist in storage
    try:
        storage.head_object(bucket=pipeline_params.storage_bucket, key=map_index_key)
        storage.head_object(bucket=pipeline_params.storage_bucket, key=filtered_map_key)
        # If they exist, return the keys and skip computing this chunk
        stats.stop_timer("function")
        return (mapper_id, map_index_key, filtered_map_key), stats
    except StorageNoSuchKeyError:
        # If any output is missing, proceed
        pass

    # Make temp dir and ch into it, save pwd to restore it later
    tmp_dir = tempfile.mkdtemp()
    pwd = os.getcwd()
    os.chdir(tmp_dir)
    print("Working directory: ", os.getcwd())
    try:
        # Get fastq chunk and store it to disk in tmp directory
        fastq_chunk_filename = f"chunk_{fastq_chunk['chunk_id']}.fastq"
        with stats.timeit("fetch_fastq_chunk"):
            fetch_fastq_chunk(pipeline_params, fastq_chunk, fastq_chunk_filename, storage)
        stats.set_value("fastq_chunk_size", os.path.getsize(fastq_chunk_filename))

        # Fetch gem file and store it to disk in tmp directory
        gem_index_filename = os.path.join(f"chunk_{fasta_chunk['chunk_id']}.gem")
        with stats.timeit("fetch_gem_chunk"):
            fetch_gem_chunk(pipeline_params, fasta_chunk, gem_index_filename, storage)
        stats.set_value("gem_chunk_size", os.path.getsize(gem_index_filename))

        # GENERATE ALIGNMENT AND ALIGNMENT INDEX (FASTQ TO MAP)
        # TODO refactor bash script
        # TODO support implement paired-end, replace not-used with 2nd fastq chunk
        # TODO use proper tmp directory instead of uuid base name

        # s3 or SRA works the same for the script /function/bin/map_index_and_filter_map_file_cmd_awsruntime.sh on single end sequences.
        cmd = [
            "/function/bin/map_index_and_filter_map_file_cmd_awsruntime.sh",
            gem_index_filename,
            fastq_chunk_filename,
            "not-used",
            pipeline_params.sra_accession,
            "s3",
            "single-end",
            str(pipeline_params.gem_mapper_threads or multiprocessing.cpu_count()),
        ]
        print(" ".join(cmd))
        with stats.timeit("map_index_and_filter"):
            out = sp.run(cmd, capture_output=True)
        print(out.stdout.decode("utf-8"))
        print(os.listdir())

        # Reorganize file names
        map_index_filename = os.path.join(tmp_dir, pipeline_params.sra_accession + "_map.index.txt")
        shutil.move(pipeline_params.sra_accession + "_map.index.txt", map_index_filename)
        filtered_map_filename = os.path.join(
            tmp_dir,
            pipeline_params.sra_accession + "_" + str(mapper_id) + "_filt_wline_no.map",
        )
        shutil.move(pipeline_params.sra_accession + "_filt_wline_no.map", filtered_map_filename)
        stats.set_value("map_index_size", os.path.getsize(map_index_filename))
        stats.set_value("filtered_map_size", os.path.getsize(filtered_map_filename))

        # Compress outputs
        zipped_map_index_filename = map_index_filename + ".bz2"
        with zipfile.ZipFile(
            zipped_map_index_filename,
            "w",
            compression=zipfile.ZIP_BZIP2,
            compresslevel=9,
        ) as zf:
            zf.write(map_index_filename, arcname=PurePosixPath(map_index_filename).name)
        stats.set_value("zipped_map_index_size", os.path.getsize(zipped_map_index_filename))

        zipped_filtered_map_filename = filtered_map_filename + ".bz2"
        with zipfile.ZipFile(
            zipped_filtered_map_filename,
            "w",
            compression=zipfile.ZIP_BZIP2,
            compresslevel=9,
        ) as zf:
            zf.write(filtered_map_filename, arcname=PurePosixPath(filtered_map_filename).name)
        stats.set_value("zipped_filtered_map_size", os.path.getsize(zipped_filtered_map_filename))

        # Copy result files to storage for index correction
        with stats.timeit("upload_map_index"):
            storage.upload_file(
                file_name=zipped_map_index_filename, bucket=pipeline_params.storage_bucket, key=map_index_key
            )

        with stats.timeit("upload_filtered_map"):
            storage.upload_file(
                file_name=zipped_filtered_map_filename, bucket=pipeline_params.storage_bucket, key=filtered_map_key
            )

        stats.stop_timer("function")
        return (mapper_id, map_index_key, filtered_map_key), stats
    finally:
        print("Cleaning up")
        os.chdir(pwd)
        shutil.rmtree(tmp_dir)


def index_correction(
    pipeline_params: PipelineParameters,
    run_id: str,
    mapper_id: str,
    map_index_keys: Tuple[str],
    storage: Storage,
):
    """
    Lithops callee function
    Corrects the index after the first map iteration.
    All the set files must have the prefix "map_index_files/".
    Corrected indices will be stored with the prefix "corrected_index/".
    """
    stats = Stats()
    stats.set_value("mapper_id", mapper_id)
    stats.start_timer("function")

    mapper_storage_tmp_prefix = partial(get_storage_tmp_prefix, run_id, "index_correction", mapper_id)

    pwd = os.getcwd()

    # TODO replace with a proper file name (maybe indicating fastq chunk id)
    output_file = "merged_filtered_index.txt"
    zipped_output_file = output_file + ".bz2"

    corrected_index_key = mapper_storage_tmp_prefix(zipped_output_file)

    # Check if output files already exist in storage
    try:
        storage.head_object(bucket=pipeline_params.storage_bucket, key=corrected_index_key)
        stats.stop_timer("function")
        # If they exist, return the keys and skip computing this chunk
        return (mapper_id, corrected_index_key), stats
    except StorageNoSuchKeyError:
        # If the output is missing, proceed
        pass

    # Download all map files for this fastq chunk
    input_temp_dir = tempfile.mkdtemp()
    output_temp_dir = tempfile.mkdtemp()
    try:
        for i, map_index_key in enumerate(map_index_keys):
            map_index_stat = Stats()
            map_index_stat.set_value("map_index_key", map_index_key)

            local_compressed_map_path = os.path.join(input_temp_dir, f"map_{i}.map.bz2")

            with map_index_stat.timeit(f"download_map_index"):
                storage.download_file(
                    bucket=pipeline_params.storage_bucket, key=map_index_key, file_name=local_compressed_map_path
                )
            map_index_stat.set_value("map_index_size", os.path.getsize(local_compressed_map_path))

            with zipfile.ZipFile(local_compressed_map_path, "r", compression=zipfile.ZIP_BZIP2, compresslevel=9) as zf:
                with map_index_stat.timeit(f"extract_map_index"):
                    zf.extractall(input_temp_dir)

            # TODO set proper map index file name
            os.rename(
                os.path.join(input_temp_dir, f"{pipeline_params.sra_accession}_map.index.txt"),
                os.path.join(input_temp_dir, f"{i}_map.index.txt"),
            )
            os.remove(local_compressed_map_path)

            stats.set_value(map_index_key, map_index_stat.dump_dict())
        print(os.listdir(input_temp_dir))

        # TODO chdir required or binary_reducer.sh
        os.chdir(output_temp_dir)

        # Execute correction scripts
        intermediate_file = f"{mapper_id}.intermediate.txt"
        # TODO binary_reducer.sh could be more efficiently implemented using python and consuming files from storage as a generator
        cmd = f"/function/bin/binary_reducer.sh /function/bin/merge_gem_alignment_metrics.sh 4 {input_temp_dir}/* > {intermediate_file}"
        print(cmd)
        with stats.timeit("merge_gem_alignment_metrics"):
            proc = sp.run(cmd, shell=True, universal_newlines=True, capture_output=True)
        print(proc.stdout)
        print(proc.stderr)
        proc.check_returncode()

        # timestamps.store_size_data("filter_merged", time())
        cmd = f"/function/bin/filter_merged_index.sh {intermediate_file} {output_file}"
        print(cmd)
        with stats.timeit("filter_merged_index"):
            proc = sp.run(cmd, shell=True, check=True, universal_newlines=True)
        print(proc.stdout)
        print(proc.stderr)
        proc.check_returncode()

        print(os.listdir())

        stats.set_value("output_file_size", os.path.getsize(output_file))

        # Compress output
        with zipfile.ZipFile(zipped_output_file, "w", compression=zipfile.ZIP_BZIP2, compresslevel=9) as zf:
            with stats.timeit("compress_output"):
                zf.write(output_file, arcname=output_file)
        stats.set_value("zipped_output_file_size", os.path.getsize(zipped_output_file))

        # Upload corrected index to storage
        with stats.timeit("upload_corrected_index"):
            storage.upload_file(
                bucket=pipeline_params.storage_bucket, key=corrected_index_key, file_name=zipped_output_file
            )

        os.chdir(pwd)
        stats.stop_timer("function")
        return (mapper_id, corrected_index_key), stats
    finally:
        os.chdir(pwd)
        force_delete_local_path(input_temp_dir)
        force_delete_local_path(output_temp_dir)


def filtered_index_to_mpileup(
    pipeline_params: PipelineParameters,
    run_id: str,
    mapper_id: str,
    fasta_chunk: dict,
    filtered_map_key: str,
    corrected_index_key: str,
    storage: Storage,
):
    stats = Stats()
    stats.set_value("mapper_id", mapper_id)
    stats.start_timer("function")

    temp_dir = tempfile.mkdtemp()
    mapper_storage_tmp_prefix = partial(get_storage_tmp_prefix, run_id, "filtered_index_to_mpileup", mapper_id)
    pwd = os.getcwd()
    # TODO get base name from params
    corrected_map_file = f"{pipeline_params.sra_accession}_{mapper_id}_filt_wline_no_corrected.map"
    mpileup_file = corrected_map_file + ".mpileup"
    mpileup_key = mapper_storage_tmp_prefix(mpileup_file)

    # Check if output file already exists in storage
    try:
        storage.head_object(bucket=pipeline_params.storage_bucket, key=mpileup_key)
        # If it exist, return the keys and skip computing this chunk
        stats.stop_timer("function")
        return (mapper_id, mpileup_key), stats
    except StorageNoSuchKeyError:
        # If the output is missing, proceed
        pass

    try:
        os.chdir(temp_dir)

        # Recover fasta chunk
        fasta_chunk_filename = f"chunk_{fasta_chunk['chunk_id']}.fasta"
        with stats.timeit("fetch_fasta_chunk"):
            fetch_fasta_chunk(fasta_chunk, fasta_chunk_filename, storage, pipeline_params.fasta_path)
        stats.set_value("fasta_chunk_size", os.path.getsize(fasta_chunk_filename))

        # Recover filtered map file
        bz2_filt_map_filename = pathlib.PurePosixPath(filtered_map_key).name
        with stats.timeit("download_filtered_map"):
            storage.download_file(
                bucket=pipeline_params.storage_bucket, key=filtered_map_key, file_name=bz2_filt_map_filename
            )
        stats.set_value("bz2_filt_map_size", os.path.getsize(bz2_filt_map_filename))

        with zipfile.ZipFile(bz2_filt_map_filename, "r", compression=zipfile.ZIP_BZIP2, compresslevel=9) as zf:
            with stats.timeit("extract_filtered_map"):
                zf.extractall()
        os.remove(bz2_filt_map_filename)

        # Get corrected map index for this fastq chunk
        try:
            bz2_corrected_index_filename = pathlib.PurePosixPath(corrected_index_key).name
        except:
            raise ValueError(corrected_index_key)

        with stats.timeit("download_corrected_index"):
            storage.download_file(
                bucket=pipeline_params.storage_bucket, key=corrected_index_key, file_name=bz2_corrected_index_filename
            )
        stats.set_value("bz2_corrected_index_size", os.path.getsize(bz2_corrected_index_filename))
        with zipfile.ZipFile(bz2_corrected_index_filename, "r", compression=zipfile.ZIP_BZIP2, compresslevel=9) as zf:
            with stats.timeit("extract_corrected_index"):
                zf.extractall()
        os.remove(bz2_corrected_index_filename)

        print(os.listdir(temp_dir))

        # Filter aligments with corrected map file
        # timestamps.store_size_data("map_file_index_correction", time())
        filt_map_filename = f"{pipeline_params.sra_accession}_{mapper_id}_filt_wline_no.map"
        # TODO replace with a proper file name (maybe indicating fastq chunk id)
        corrected_index_filename = "merged_filtered_index.txt"

        cmd = [
            "/function/bin/map_file_index_correction.sh",
            corrected_index_filename,
            filt_map_filename,
            str(pipeline_params.tolerance),
        ]
        print(" ".join(cmd))
        with stats.timeit("map_file_index_correction"):
            proc = sp.run(cmd, capture_output=True)  # change to _v3.sh and runtime 20
        print(proc.stdout.decode("utf-8"))
        print(proc.stderr.decode("utf-8"))

        print(os.listdir(temp_dir))

        # Generate mpileup
        # timestamps.store_size_data("gempileup_run", time())
        cmd = [
            "/function/bin/gempileup_run.sh",
            corrected_map_file,
            fasta_chunk_filename,
        ]
        print(" ".join(cmd))
        with stats.timeit("gempileup_run"):
            proc = sp.run(cmd, capture_output=True)
        print(proc.stdout.decode("utf-8"))
        print(proc.stderr.decode("utf-8"))

        # Store output to storage
        stats.set_value("mpileup_size", os.path.getsize(mpileup_file))
        with stats.timeit("upload_mpileup"):
            storage.upload_file(bucket=pipeline_params.storage_bucket, key=mpileup_key, file_name=mpileup_file)

        stats.stop_timer("function")
        return (mapper_id, mpileup_key), stats
    finally:
        os.chdir(pwd)
        force_delete_local_path(temp_dir)


def mpileup_conversion(
    self,
    mpileup_file: str,
    fasta_chunk: dict,
    fastq_chunk: str,
    exec_param: str,
    storage: Storage,
) -> Tuple[str]:
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
    with open(mpileup_file, "r") as f:
        rows = f.read().splitlines()
        content = [row.split("\t") for row in rows]
        content.pop(-1)
        del rows

    # Convert mpileup to Pandas dataframe
    df = pd.DataFrame(data=content)
    df.columns = df.columns.astype(str)
    df["1"] = df["1"].astype(int64)

    # Remove disallowed characters
    fasta_key = self.fasta_chunks_prefix
    disallowed_characters = "._-!/·¡"
    for character in disallowed_characters:
        fasta_key = fasta_key.replace(character, "")

    # Create intermediate key
    fasta_chunk = str(fasta_chunk["id"])
    max_index = df.iloc[-1]["1"]
    intermediate_key = (
        self.args.file_format
        + "/"
        + exec_param
        + "/"
        + fasta_key
        + "_"
        + fasta_chunk
        + "-"
        + fastq_chunk[0]
        + "_chunk"
        + str(fastq_chunk[1]["number"])
        + "_"
        + str(max_index)
    )

    range_index = []
    x = 0
    while x < max_index:
        if x < max_index:
            range_index.append(x)
        x = x + 100000

    x = x + (max_index - x)
    range_index.append(x)

    df3 = df.groupby(pd.cut(df["1"], range_index)).count()
    content = ""
    for i in range(len(df3)):
        content = content + str(range_index[i + 1]) + ":" + str(df3.iloc[i, 0]) + "\n"

    # Upload .txt file to storage
    storage.put_object(bucket=self.args.storage_bucket, key=intermediate_key + ".txt", body=content)

    # Write the mpileup file to the tmp directory
    if self.args.file_format == "csv":
        df.to_csv(mpileup_file + ".csv", index=False, header=False)
        with open(mpileup_file + ".csv", "rb") as f:
            storage.put_object(bucket=self.args.storage_bucket, key=intermediate_key + ".csv", body=f)

    elif self.args.file_format == "parquet":
        df.to_parquet(mpileup_file + ".parquet")
        with open(mpileup_file + ".parquet", "rb") as f:
            storage.put_object(
                bucket=self.args.storage_bucket,
                key=intermediate_key + ".parquet",
                body=f,
            )

    return [intermediate_key + "." + self.args.file_format, intermediate_key + ".txt"]
