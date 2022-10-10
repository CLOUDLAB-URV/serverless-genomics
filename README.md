# SERVERLESS VARIANT CALLER PIPELINE - USING LITHOPS + GEM-MAPPER + SiNPle

## INSTALLATION REQUIREMENTS 

### 1. Install local dependencies (where the script is executed):

- python 3.8 and `pip install lithops\[aws\]`
- redis (redis-cli)
- edirect (if using SRA option): `sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"`. Important to have the new script added binaries to the PATH variable of the shell used to launch the variant caller.


### 2. Build and upload the runtime

`cd dockerfile && lithops runtime build -f Dockerfile -b aws_lambda lumimar/hutton-genomics-v03:13`


## RUNNING THE PIPELINE

This is a command line example in which all configurable parameters are made explicit, although many have default values that can be used so that only few parameters are actually required.

```
python varcall_lithops_demo_v5.py -fq ERR9729866 -fa hg19.fa -cl aws -b cloudbutton-variant-caller-input -fb ayman-lithops-meta-cloudbutton-hutton -ds SRA -nfq 2000000 -nfa 100000000 -ofa 300 -rl 152 -t 0 -ff csv -s3w False -rt lumimar/hutton-genomics-v03:18 -rtm 4096 -rtr 4096 -bs 75% -ftm 900 -ftr 900 -sk False -ip 54.146.89.181 -id i-0cee52f66655d990b -rg us-east-1

```

A wrapper script, `run_variant_caller_v2.sh`, allows to run the variant caller with these additional features:
- inputting command line arguments from a table `varcall_args.tsv`
- sorting live function log by function number for readability
- extraction of summary stats for the run, to quickly identify where the run might be failing
- summary stats linked to log file and to command line arguments
- generation of plots with time of various pipeline stages and size of intermediate files generated
- specifying command and run id, and how many repeats of the script to run, if testing reproducibility. 

Here is a command line example:
```
bash run_variant_caller_v2.sh varcall_args.tsv 1 1 True
```
The four command lines arguments are: 1. variant caller argument table, 2. run id (if same settings run multiple times), 3. number of iterations (if script run multiple times), 4. whether to run the variant caller (True) or just process the log file (False). The latter is useful in the event of an aborted run, as it processes the existing incomplete log and generates plots and info to quickly identify where a problem arose.


Together with the sorted log and a series of tables (all saved to the `varcall_out` subfolder), the script also adds lines to the summary file `varcall_results_summary.tsv`, providing basic stats about the number of functions that "made it" through the various stages of the pipeline. Its output also includes the initial commandline, and is also present at the end of the log file.

The args `varcall_args.tsv` has two header lines, with short and long argument option names (presented here in column format, first two columns in table below)

short option name | long option name | explanation | example
------------------|------------------|-------------|--------
run_n | #NA | NA | #1
fq | fq_seq_name | fastq sequence name | SRR15068323
fq1 | fastq1 | fastq file name 1 | NA
fq2 | fastq2 | fastq file name 2 | NA
fa | fasta | fasta file name | hg19.fa
cl | cloud_adr | cloud provider | aws
b | bucket | bucket | cloudbutton-variant-caller-input
fb | fbucket | bucket for fasta file | ayman-lithops-meta-cloudbutton-hutton
ds | data_source | SRA or s3  | SRA
nfq | fastq_read_n | number of reads per fastq chunk | 800000
nfa | fasta_char_n | number of characters per fasta chunk | 100000000
ofa | fasta_char_overlap | overlap between fasta chunks | 300
rl | read_length | read length (to calculate approx fastq size) | 152
t | tolerance | number of additional strata to include in alignment output | 0
ff | file_format | mpileup file format conversion choice (csv or parquet) | csv
itn | iterdata_n | Number of functions to launch (all if empty) | 5
s3w | temp_to_s3 | saving temporary files to s3 for debugging | FALSE
rt | runtime_id | lambda function runtime id | lumimar/hutton-genomics-v03:18
rtm | runtime_mem | memory associated with map function  | 4096
rtr | runtime_memr | memory associated with reduce function | 4096
bs | buffer_size | size of buffer in reduce function | 75%
ftm | func_timeout_map | timeout for map functions | 900
ftr | func_timeout_reduce | timeout for reduce functions | 900
sk | skip_map | skip map function if debugging reduce phase | FALSE
ip | ec2ip | ec2 IP for redis | 54.146.89.181
id | ec2id | ec2 id for redis | i-0cee52f66655d990b
rg | ec2region | ec2 region | us-east-1

```