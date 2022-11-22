from dataclasses import dataclass


def validate_parameters(params):
    if 'fasta_bucket' not in params:
        params['fasta_bucket'] = params['bucket']

    params['fastq_chunk_size'] = 4 * int(params['fastq_read_n'])

    if 'fastq2' not in params:
        params['seq_type'] = "single-end"
    else:
        params['seq_type'] = "paired-end"

    if params['checkpoints'] in ('True', 'true', 'Yes', 'yes', 'y', '1'):
        params['checkpoints'] = True
    else:
        params['checkpoints'] = False

    params['execution_name'] = f"{params['fq_seqname']}_" \
                               f"{params['fasta_file']}_" \
                               f"{params['fastq_read_n']}_" \
                               f"{params['fasta_workers']}"

    return PipelineParameters(**params)


@dataclass(frozen=True)
class PipelineParameters:
    # Mandatory Parameters
    fasta_file: str
    bucket: str
    fasta_bucket: str
    fasta_workers: int
    fastq_chunk_size: int
    seq_type: str
    execution_name: str

    # Optional Parameters
    # FastQ names
    fq_seqname: str = None
    fastq_file: str = ""

    # From input, determine whether it is paired- or single-end sequencing
    fastq_file2: str = ""

    # Cloud Storage Settings
    cloud_adr: str = "aws"

    # Fastq data source (SRA)
    datasource: str = "s3"

    # File Splitting Parameters
    fastq_read_n: int = None

    # Pipeline-Specific Parameters
    tolerance: int = 0
    file_format: str = "parquet"

    # Run Settings
    iterdata_n: ... = None
    function_n: ... = "All"
    concur_fun: int = 10000
    temp_to_s3: bool = False
    runtime_id: str = 'lumimar/hutton-genomics-v03:18'
    runtime_mem: int = 1024
    runtime_mem_r: int = 4096
    runtime_storage: int = 4000
    buffer_size: int = "75%"
    func_timeout_map: int = 2400
    func_timeout_reduce: int = 2400
    skip_map: bool = False
    lb_method: str = "select"
    checkpoints: bool = False

    # Debug Parameters
    gem_test: bool = False
    pre_processing_only: bool = False
    debug: bool = True

    # S3 prefixes (cloud folders)
    fasta_folder: str = "fasta/"
    fastq_folder: str = "fastqgz/"
    fasta_folder_index: str = "fasta-indexes/"
    idx_folder: str = "fastq-indexes/"
    out_folder: str = "outputs/"
    s3_temp_folder: str = "temp_outputs/"
