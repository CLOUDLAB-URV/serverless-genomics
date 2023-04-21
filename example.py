import logging

from serverlessgenomics import VariantCallingPipeline

BUCKET = "serverless-genomics"

if __name__ == '__main__':
    pipeline = VariantCallingPipeline(
        fasta_path=f's3://{BUCKET}/fasta/TriTrypDB-9.0_TbruceiTREU927_Genome_MbChr.fasta',
        fasta_chunks=5,
        fastq_path=f's3://{BUCKET}/fastq/SRR935389.fastq.gz',
        fastq_chunks=20,
        log_level=logging.DEBUG,
        storage_bucket=BUCKET,
        log_stats=True,
        debug=True,
        gem_mapper_threads=1
    )
    # pipeline.clean_all()
    # pipeline.preprocess()
    pipeline.run_pipeline()
