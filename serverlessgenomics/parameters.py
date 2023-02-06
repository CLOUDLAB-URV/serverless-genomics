from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import Optional

import lithops
import uuid

from .cachedlithops import CachedLithopsInvoker
from .utils import S3Path


@dataclass(frozen=True)
class PipelineRun:
    """
    Dataclass to store a pipeline's execution state
    """
    # Storage path for fasta input file
    fasta_path: S3Path
    # Number of chunks to split fasta input file into
    fasta_chunks: int
    # Number of chunks to split fastq input file into
    fastq_chunks: int
    # Storage path for fastq input file
    fastq_path: Optional[S3Path] = None
    # SRA identifier for fastq read input
    fastq_sra: Optional[str] = None
    # Execution run UUID
    run_id: str = str(uuid.uuid4())
    # TODO what is tolerance? (ask Lucio)
    tolerance: int = 0
    # fastq chunks to be processed
    fastq_chunk_range: range = None
    # FASTQ base name
    base_name: str = "SRRXXXXXX"

    # Lithops settings
    max_workers: int = 1000
    runtime_image: str = None
    runtime_mem: int = 1024
    runtime_timeout: int = 2400
    func_timeout_reduce: int = 2400
    lb_method: str = 'select'
    checkpoints: bool = False
    log_level: int = logging.INFO
    log_stats: bool = False
    log_stats_name: str = 'logs_stats'

    # Bucket name with write permissions to store preprocessed, intermediate and output data
    storage_bucket: str = 'serverless-genomics'
    # Prefix for fastqgz generated index keys
    fastqgz_idx_prefix: str = 'fastqgz-indexes/'
    # Prefix for faidx generated index keys
    faidx_prefix: str = 'faidx-indexes/'
    # Prefix for generated reference genome indexes
    genome_index_prefix: str = 'fasta-indexes/'
    # Prefix for output data keys
    output_prefix: str = 'output/'
    # Prefix for temporal data keys
    tmp_prefix: str = 'tmp/'

    @property
    def fastqgz_idx_keys(self):
        """
        Returns a tuple for fastqgz index keys in storage as (index file key, tab data key)
        """
        return (os.path.join(self.fastqgz_idx_prefix, self.fastq_path.key + '.idx'),
                os.path.join(self.fastqgz_idx_prefix, self.fastq_path.key + '.tab'))

    @property
    def faidx_key(self):
        return os.path.join(self.faidx_prefix, self.fasta_path.key + '.fai')


@dataclass(frozen=True)
class Lithops:
    storage: lithops.Storage
    invoker: CachedLithopsInvoker


def validate_parameters(params) -> PipelineRun:
    if 'fasta_path' not in params:
        raise KeyError()
    if 'fasta_chunks' not in params:
        raise KeyError()

    try:
        params['run_id'] = params['override_id']
        del params['override_id']
    except:
        pass

    params['fasta_path'] = S3Path.from_uri(params['fasta_path'])
    if 'fastq_path' in params:
        params['fastq_path'] = S3Path.from_uri(params['fastq_path'])

    fastq_file = str(params['fastq_path']).split("/")[-1]
    if(fastq_file[0:3] == "SRR"):
        params["base_name"] = fastq_file.split(".")[0]

    return PipelineRun(**params)
