from __future__ import annotations

import multiprocessing
import os
import subprocess
import tempfile
from time import time
from typing import TYPE_CHECKING


from .. import fetch_fasta_chunk
from ...stats import Stats
from ...utils import force_delete_local_path

if TYPE_CHECKING:
    from ...pipeline import PipelineParameters
    from typing import List
    from lithops import Storage


def get_gem_chunk_storage_key(pipeline_params: PipelineParameters, fasta_chunk_id: int) -> str:
    return os.path.join(
        pipeline_params.gem_index_prefix,
        pipeline_params.fasta_path.key,
        f"{pipeline_params.fasta_chunks}-chunks",
        "chunk" + str(fasta_chunk_id).zfill(4) + ".gem",
    )


def get_gem_chunk_storage_prefix(pipeline_params: PipelineParameters) -> str:
    return os.path.join(
        pipeline_params.gem_index_prefix, pipeline_params.fasta_path.key, f"{pipeline_params.fasta_chunks}-chunks"
    )


def gem_indexer(pipeline_params: PipelineParameters, fasta_chunk_id: int, fasta_chunk: dict, storage: Storage):
    # Stats
    stat, timestamps, data_size = Stats(), Stats(), Stats()
    stat.timer_start(fasta_chunk_id)
    timestamps.store_size_data("start", time())

    # Initialize names
    gem_index_filename = os.path.join(f"chunk{str(fasta_chunk_id).zfill(4)}.gem")

    # Check if gem file already exists
    # try:
    #     storage.head_object(bucket=pipeline_params.storage_bucket, key=f"gem/{gem_index_filename}")
    #     stat.timer_stop(fasta_chunk_id)
    #     return stat.get_stats()
    # except StorageNoSuchKeyError:
    #     pass

    # Make temp dir and ch into it, save cwd to restore it later
    tmp_dir = tempfile.mkdtemp()
    cwd = os.getcwd()
    os.chdir(tmp_dir)

    try:
        # Get fasta chunk and store it to disk in tmp directory
        timestamps.store_size_data("download_fasta", time())

        fasta_chunk_filename = f"chunk_{str(fasta_chunk['chunk_id']).zfill(4)}.fasta"
        fetch_fasta_chunk(fasta_chunk, fasta_chunk_filename, storage, pipeline_params.fasta_path)

        data_size.store_size_data(fasta_chunk_filename, os.path.getsize(fasta_chunk_filename) / (1024 * 1024))

        # gem-indexer already appends .gem to output file, so we delete it from the process call arguments
        timestamps.store_size_data("gem_indexer", time())
        cmd = [
            "gem-indexer",
            "--input",
            fasta_chunk_filename,
            "--threads",
            str(multiprocessing.cpu_count()),
            "-o",
            gem_index_filename.replace(".gem", ""),
        ]
        print(" ".join(cmd))
        out = subprocess.run(cmd, capture_output=True)
        print(out.stderr.decode("utf-8"))

        # TODO apparently 1 is return code for success (why)
        assert out.returncode == 1

        # Upload the gem file to storage
        timestamps.store_size_data("upload_gem", time())
        gem_index_key = get_gem_chunk_storage_key(pipeline_params, fasta_chunk_id)
        storage.upload_file(file_name=gem_index_filename, bucket=pipeline_params.storage_bucket, key=gem_index_key)

        # Finish stats
        data_size.store_size_data(gem_index_filename, os.path.getsize(gem_index_filename) / (1024 * 1024))
        stat.timer_stop(fasta_chunk_id)
        timestamps.store_size_data("end", time())

        # Store dict
        stat.store_dictio(timestamps.get_stats(), "timestamps", fasta_chunk_id)
        stat.store_dictio(data_size.get_stats(), "data_sizes", fasta_chunk_id)
        return stat.get_stats()
    finally:
        os.chdir(cwd)
        force_delete_local_path(tmp_dir)
