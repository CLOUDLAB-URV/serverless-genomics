# Serverless Variant Calling

## Datasets

1. Trypanosome 

- Reference genome (FASTA, 35 MB): (Link)[https://tritrypdb.org/tritrypdb/app/downloads/Current_Release/TbruceiTREU927/fasta/data/]
- Sequence Reads:
    * SRR6052133* (FASTQGZip, 668 MB): (Link)[https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR6052133&display=download]

2. Human

- Reference genome (FASTA GZip, 905 MB): (Link)[http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/]
- Sequence Reads:
    * SRR15068323 (FASTQGZip, 1.2 GB): (Link)[https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR15068323&display=data-access]
    * ERR9856489 (FASTQGZip, 12.1 GB): (Link)[https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR9856489&display=data-access]

3. Bos taurus

- Reference genome (FASTA GZip, 781 MB): (Link)[https://www.ensembl.org/Bos_taurus/Info/Index]
- Sequence Reads:
    * SRR934415* (FASTQGZip, 16.5 GB): (Link)[https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR934415&display=data-access]


\*The pipeline currently support only single-end sequences. Use only the first sequence for paired-end reads.

## Setup

1. Setup Lithops for AWS backend.

2. Build the runtime in the `dockerfile` directory :

```bash
$ lithops runtime build -f lambda.Dockerfile serverless-genomics:1
```

3. Configure Lithops to use the built runtime (e.g. `serverless-genomics:1`). The required runtime memory, ephemeral disk size and timeout will depend on the number of partitions and 

4. Create an S3 bucket with the required input datasets.

5. Use `VariantCallingPipeline` to execute the pipeline. Provide the necessary parameters. You can see a complete list in file `serverlessgenomics/pipeline.py`. The required parameters are:
    - `run_id`: ID of a specific run. It can be reused for failed runs. You must change the ID if different input data are used.
    - `fasta_path`: Path of the input sequence read (`s3://...`).
    - `fasta_chunks`: Number of FASTA partitions.
    - `fastq_path`: Path of the input reference genome (`s3://...`).
    - `fastq_chunks`: Number of FASTQ partitions.
    - `storage_bucket`: Temporary data S3 bucket name. It must exist.

6. Call `VariantCallingPipeline.preprocess()` for the pre-processing stage, `VariantCallingPipeline.alignment()` to run the alignment phase, `VariantCallingPipeline.reduce()` for the alignment phase or `VariantCallingPipeline.run_pipeline()` for a complete execution.