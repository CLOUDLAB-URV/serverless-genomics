from dataclasses import dataclass
from re import I

@dataclass
class Arguments:
    # FastQ names
    fq_seqname: str
    fastq_file: str
    
    # From input, determine whether it is paired- or single-end sequencing
    fastq_file2: str
    seq_type: str
    fasta_file: str
    
    # Cloud Storage Settings
    cloud_adr: str
    bucket: str
    fasta_bucket: str
    
    # Fastq data source (SRA)
    datasource: str
    
    # File Splitting Parameters
    fastq_read_n: int
    fastq_chunk_size: int
    fasta_chunk_size: int
    fasta_char_overlap: int
    
    # Pipeline-Specific Parameters
    tolerance: int
    file_format: str
    
    # Run Settings
    iterdata_n: ...
    function_n: ...
    concur_fun: int
    temp_to_s3: bool
    runtime_id: str
    runtime_mem: int
    runtime_mem_r: int
    runtime_storage: int
    buffer_size: int
    func_timeout_map: int
    func_timeout_reduce: int
    skip_map: bool
    lb_method: str
    
    # DEBUGGING SETTINGS
    gem_test: bool
    pre_processing_only: bool
    debug: bool
    
    # S3 prefixes (cloud folders)
    fasta_folder: str
    fastq_folder: str
    split_fasta_folder: str
    idx_folder: str
    out_folder: str
    s3_temp_folder: str