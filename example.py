import logging

from serverlessgenomics import VariantCallingPipeline

if __name__ == '__main__':
    pipeline = VariantCallingPipeline(
        fasta_path='s3://agabriel-data/fasta/TriTrypDB-9.0_TbruceiTREU927_Genome_MbChr.fasta',
        fasta_chunks=5,
        fastq_path='s3://agabriel-data/fastq/SRR935389.fastq.gz',
        fastq_chunks=20,
        log_level=logging.DEBUG,
        storage_bucket='agabriel-data',
        #fastq_chunk_range=range(10,13),
        override_id='616d81c4-192c-4775-92f0-41bcf3653d0a'
    )
    # pipeline = VariantCallingPipeline.restore_run('616d81c4-192c-4775-92f0-41bcf3653d0a')
    # pipeline.preprocess()
    pipeline.run_pipeline()
