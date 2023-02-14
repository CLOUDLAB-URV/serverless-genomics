import logging

from serverlessgenomics import VariantCallingPipeline

if __name__ == '__main__':
    pipeline = VariantCallingPipeline(
        fasta_path='s3://agabriel-data/fasta/TriTrypDB-9.0_TbruceiTREU927_Genome_MbChr.fasta',
        fasta_chunks=5,
        fastq_path='s3://agabriel-data/fastq/SRR6052133_1.fastq.gz',
        fastq_chunks=20,
        log_level=logging.DEBUG,
        storage_bucket='agabriel-data',
        override_id='8b3cead0-fd3f-4b1b-812e-618b6dbc4a28'
    )
    pipeline.run_pipeline()
