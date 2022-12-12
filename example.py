import logging

from serverlessgenomics import VariantCallingPipeline

if __name__ == '__main__':
    pipeline = VariantCallingPipeline(
        fasta_path='s3://aitor-serverless-genomics/fasta/TriTrypDB-9.0_TbruceiTREU927_Genome_MbChr.fasta',
        fasta_chunks=4,
        fastq_path='s3://aitor-serverless-genomics/fastq/SRR935389.fastq.gz',
        fastq_chunks=4,
        log_level=logging.DEBUG,
        storage_bucket='aitor-serverless-genomics',
        override_id='616d81c4-192c-4775-92f0-41bcf3653d0a'
    )
    # pipeline = VariantCallingPipeline.restore_run('616d81c4-192c-4775-92f0-41bcf3653d0a')
    # pipeline.preprocess()
    pipeline.run_pipeline()
