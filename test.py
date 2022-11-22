import logging

from serverlessgenomics import VariantCallingPipeline

if __name__ == '__main__':
    pipeline = VariantCallingPipeline(
        fasta_file='TriTrypDB-9.0_TbruceiTREU927_Genome_MbChr.fasta',
        bucket='genomics-pipeline',
        fq_seqname='SRR935389',
        fastq_read_n=4,
        fasta_workers=4,
        log_level=logging.DEBUG,
        fasta_folder=""
    )
    pipeline.run_pipeline()
