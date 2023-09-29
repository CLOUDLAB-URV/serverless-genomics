from __future__ import annotations

import logging
import multiprocessing
import os
import re
import subprocess
import tempfile
from time import time
from typing import TYPE_CHECKING, Set
from lithops.storage.utils import StorageNoSuchKeyError

from ..datasource import fetch_fasta_chunk
from ..datasource.sources.gem import (
    get_gem_chunk_storage_key,
    get_gem_chunk_storage_prefix,
)
from ..utils import force_delete_local_path
from ..stats import Stats

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

    # Find cached gem files for this FASTA file and chunk size
    cached_gems_keys = lithops.storage.list_keys(bucket=pipeline_params.storage_bucket, prefix=gems_prefix)
    cached_gem_chunk_ids = []

    for cached_gems_key in cached_gems_keys:
        basename = os.path.basename(cached_gems_key)
        matches = re.findall(r"\d+", basename)
        assert len(matches) == 1
        chunk_id = int(matches.pop())
        cached_gem_chunk_ids.append(chunk_id)

    cached_gem_chunk_ids = set(cached_gem_chunk_ids)
    requested_gems_ids = {fq_ch["chunk_id"] for fq_ch in fasta_chunks}

    # Compare cached gem file set and requested gem file set
    if requested_gems_ids.issubset(cached_gem_chunk_ids):
        # All requested chunks are already in storage
        logger.info('Using %d cached GEM files in storage (prefix="%s")', len(requested_gems_ids), gems_prefix)
        return cached_gem_chunk_ids, {}

    # Generate missing gem files
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
    results = lithops.invoker.map(gem_indexer, iterdata)
    gem_keys, stats = zip(*results)

    return gem_keys, stats


def generate_gem_indexer_iterdata(pipeline_params: PipelineParameters, fasta_chunks: List[dict]) -> List[dict]:
    iterdata = []

    for fa_ch in fasta_chunks:
        params = {
            "pipeline_params": pipeline_params,
            "fasta_chunk_id": fa_ch["chunk_id"],
            "fasta_chunk": fa_ch,
        }
        iterdata.append(params)

    return iterdata


def gem_indexer(pipeline_params: PipelineParameters, fasta_chunk_id: int, fasta_chunk: dict, storage: Storage):
    stats = Stats()
    stats.set_value("fasta_chunk_id", fasta_chunk_id)

    # Initialize names
    gem_index_filename = os.path.join(f"chunk{str(fasta_chunk_id).zfill(4)}.gem")
    gem_index_key = get_gem_chunk_storage_key(pipeline_params, fasta_chunk_id)

    # Check if gem file already exists
    try:
        storage.head_object(bucket=pipeline_params.storage_bucket, key=gem_index_key)
        return gem_index_key, stats
    except StorageNoSuchKeyError:
        pass

    # Make temp dir and ch into it, save cwd to restore it later
    tmp_dir = tempfile.mkdtemp()
    cwd = os.getcwd()
    os.chdir(tmp_dir)

    try:
        # Get fasta chunk and store it to disk in tmp directory
        with stats.timeit("fetch_fasta_chunk"):
            fasta_chunk_filename = f"chunk_{str(fasta_chunk['chunk_id']).zfill(4)}.fasta"
            fetch_fasta_chunk(fasta_chunk, fasta_chunk_filename, storage, pipeline_params.fasta_path)
        stats.set_value("fasta_chunk_size", os.path.getsize(fasta_chunk_filename))

        # gem-indexer already appends .gem to output file, so we delete it from the process call arguments
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
        with stats.timeit("gem_indexer"):
            out = subprocess.run(cmd, capture_output=True)
        print(out.stderr.decode("utf-8"))

        # TODO apparently 1 is return code for success (why)
        assert out.returncode == 1

        # Upload the gem file to storage
        gem_index_key = get_gem_chunk_storage_key(pipeline_params, fasta_chunk_id)
        with stats.timeit("upload_gem_index"):
            storage.upload_file(file_name=gem_index_filename, bucket=pipeline_params.storage_bucket, key=gem_index_key)
        stats.set_value("gem_index_size", os.path.getsize(gem_index_filename))

        return gem_index_key, stats
    finally:
        os.chdir(cwd)
        force_delete_local_path(tmp_dir)
