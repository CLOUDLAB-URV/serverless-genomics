from __future__ import annotations

from lithops import Storage
import subprocess as sp
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from serverlessgenomics.parameters import PipelineRun
import pathlib

from .fasta_partitioner_index import fasta_partitioner_caller
import re


def prepare_fasta(args: PipelineRun):
    """
    Try to find fasta chunks in storage. If they are not found proceed to partition the orignal fasta file.
    """
    storage = Storage()
    fasta_index = args.fasta_indexes_prefix + pathlib.Path(args.fasta_path).stem + '.fai'
    try:
        print("Searching " + args.bucket + "/" + fasta_index)
        storage.head_object(args.bucket, fasta_index)
    except:
        print("Index fasta folder empty / not found, then generating index fasta file and storing it")
        # Generate index file from the fasta source
        fasta_partitioner_caller(args.bucket, args.fasta_folder + args.fasta_path, int(args.fasta_chunks),
                                 args.fasta_indexes_prefix)
        try:
            storage.head_object(args.bucket, fasta_index)
        except:
            print("Error generating fasta file index")
    return fasta_index


def create_fasta_chunk_for_runtime(storage: Storage, bucket: str, fasta: dict, byte_range: dict, folder: str,
                                   file_name: str):
    extra_args = {'Range': f"bytes={fasta['chunk']['offset_head']}-{fasta['chunk']['offset_base']}"}
    data = list(re.finditer(r">.+\n",
                            storage.get_object(bucket=bucket, key=folder + file_name, extra_get_args=extra_args).decode(
                                'utf-8')))[0].group()
    base = storage.get_object(bucket=bucket, key=folder + file_name, extra_get_args=byte_range).decode('utf-8')

    data += base[1::] if base[0:1] == '\n' else base  # Data has already a '\n', (data == >...\n), avoid doble '\n'

    return data.encode('utf-8')
