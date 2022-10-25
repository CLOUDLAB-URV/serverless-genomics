from dataclasses import dataclass
from re import I

@dataclass
class Arguments():
    ###################################################################
    #### MANDATORY VARIABLES
    ###################################################################
    fasta_file: str
    bucket: str
    fasta_bucket: str
    fasta_chunk_size: int
    fastq_chunk_size: int
    seq_type: str
    
    ###################################################################
    #### OPTIONAL VARIABLES
    ###################################################################
    # FastQ names
    fq_seqname : str = None
    fastq_file: str = ""
    
    # From input, determine whether it is paired- or single-end sequencing
    fastq_file2: str = ""
    
    # Cloud Storage Settings
    cloud_adr: str = "aws"
    
    # Fastq data source (SRA)
    datasource: str = "s3"
    
    # File Splitting Parameters
    fastq_read_n: int = None
    fasta_char_overlap: int = 300
    
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
    
    # DEBUGGING SETTINGS
    gem_test: bool = False
    pre_processing_only: bool = False
    debug: bool = True
    
    # S3 prefixes (cloud folders)
    fasta_folder: str = "fasta/"
    fastq_folder: str = "fastqgz/"
    split_fasta_folder: str = "fasta-chunks/"
    idx_folder: str = "fastq-indexes/"
    out_folder: str = "outputs/"
    s3_temp_folder: str = "temp_outputs/"