from __future__ import annotations

import re
from typing import TYPE_CHECKING

from .datasources import FASTQSource
from .sources.fastqgz import fetch_fastq_chunk_s3_fastqgzip
from .sources.sra import fetch_fastq_chunk_sra

if TYPE_CHECKING:
    from lithops import Storage
    from ..utils import S3Path
from ..pipelineparams import PipelineParameters


def fetch_fastq_chunk(pipeline_parameters: PipelineParameters, fastq_chunk: dict, target_filename: str,
                      storage: Storage):
    assert "source" in fastq_chunk
    if fastq_chunk["source"] == FASTQSource.S3_GZIP:
        fetch_fastq_chunk_s3_fastqgzip(
            fastq_chunk,
            target_filename,
            pipeline_parameters,
            storage
        )
    elif fastq_chunk["source"] == FASTQSource.SRA:
        fetch_fastq_chunk_sra(pipeline_params.sra_accession, fastq_chunk, fastq_chunk_filename)
    else:
        raise KeyError(fastq_chunk["source"])


def fetch_fasta_chunk(fasta_chunk: dict, target_filename: str, storage: Storage, fasta_path: S3Path):
    # Get header data
    extra_args = {"Range": f"bytes={fasta_chunk['offset_head']}-{fasta_chunk['offset_base']}"}
    header_body = storage.get_object(bucket=fasta_path.bucket, key=fasta_path.key, extra_get_args=extra_args)
    chunk_body = list(re.finditer(r">.+\n", header_body.decode("utf-8")))[0].group()

    # Get chunk body and append to header
    extra_args = {"Range": f"bytes={fasta_chunk['offset_base']}-{fasta_chunk['last_byte']}"}
    base = storage.get_object(bucket=fasta_path.bucket, key=fasta_path.key, extra_get_args=extra_args).decode("utf-8")
    chunk_body += base[1::] if base[0:1] == "\n" else base

    with open(target_filename, "w") as target_file:
        target_file.writelines(chunk_body)
