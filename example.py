import logging

from serverlessgenomics import VariantCallingPipeline

if __name__ == '__main__':
    pipeline = VariantCallingPipeline(
        fasta_path='s3://serverless-genomics/fasta/TriTrypDB-9.0_TbruceiTREU927_Genome_MbChr.fasta',
        fasta_chunks=4,
        fastq_path='s3://serverless-genomics/fastq/SRR935389.fastq.gz',
        fastq_chunks=4,
        log_level=logging.DEBUG,
    )
    # pipeline.clean_all()
    pipeline.run_pipeline()
