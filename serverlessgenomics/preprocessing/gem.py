from __future__ import annotations

import logging
import multiprocessing
import os
import re
import subprocess
import tempfile
from time import time
from typing import TYPE_CHECKING, Set

from ..datasource import fetch_fasta_chunk
from ..datasource.sources.gem import get_gem_chunk_storage_key, get_gem_chunk_storage_prefix
from ..stats import Stats
from ..utils import force_delete_local_path

if TYPE_CHECKING:
    from typing import List
    from ..pipeline import PipelineParameters, Lithops
    from lithops import Storage

logger = logging.getLogger(__name__)


def prepare_gem_chunks(pipeline_params: PipelineParameters, fasta_chunks: list[dict], lithops: Lithops) -> Set[int]:
    """
    Generate GEM indexed file metadata
    """
    # Check if all GEM files exist for the input FASTA file and specified number of chunks
    gems_prefix = get_gem_chunk_storage_prefix(pipeline_params)

    cached_gems_keys = lithops.storage.list_keys(bucket=pipeline_params.storage_bucket, prefix=gems_prefix)
    cached_gem_chunk_ids = []

    for cached_gems_key in cached_gems_keys:
        basename = os.path.basename(cached_gems_key)
        print(basename)
        matches = re.findall(r"\d+", basename)
        assert len(matches) == 1
        chunk_id = int(matches.pop())
        cached_gem_chunk_ids.append(chunk_id)

    cached_gem_chunk_ids = set(cached_gem_chunk_ids)
    requested_gems_ids = set(range(pipeline_params.fasta_chunks))

    if cached_gem_chunk_ids == requested_gems_ids:
        # All chunks are already in storage
        logger.info('Using %d cached GEM files in storage (prefix="%s")', len(cached_gem_chunk_ids), gems_prefix)
        return cached_gem_chunk_ids

    if not cached_gem_chunk_ids:
        # All chunks are missing
        iterdata = generate_gem_indexer_iterdata(pipeline_params, fasta_chunks)
        logger.debug("All GEM chunks are missing (%d in total)", len(requested_gems_ids))
    else:
        # Only some chunks are missing
        missing_chunks_ids = requested_gems_ids - cached_gem_chunk_ids
        logger.debug("Some GEM chunks are missing (%s)", missing_chunks_ids.__repr__())
        missing_fasta_chunks = []
        for fa_ch in fasta_chunks:
            if fa_ch["chunk_id"] not in missing_chunks_ids:
                missing_fasta_chunks.append(fa_ch)
        iterdata = generate_gem_indexer_iterdata(pipeline_params, missing_fasta_chunks)

    logger.info("Going to index %d GEM chunks", len(iterdata))
    lithops.invoker.map(gem_indexer, iterdata)

    return requested_gems_ids


def generate_gem_indexer_iterdata(pipeline_params: PipelineParameters, fasta_chunks: List[dict]) -> List[dict]:
    iterdata = []

    for fa_ch in fasta_chunks:
        params = {"pipeline_params": pipeline_params, "fasta_chunk_id": fa_ch["chunk_id"], "fasta_chunk": fa_ch}
        iterdata.append(params)

    return iterdata


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
