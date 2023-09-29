import logging
import json

from serverlessgenomics import VariantCallingPipeline

BUCKET = "agabriel-data"

if __name__ == "__main__":
    pipeline = VariantCallingPipeline(
        run_id="small_run",
        fasta_path=f"s3://{BUCKET}/fasta/TriTrypDB-9.0_TbruceiTREU927_Genome_MbChr.fasta",
        fasta_chunks=20,
        fastq_path=f"s3://{BUCKET}/fastq/SRR6052133_1.fastq.gz",
        fastq_chunks=10,
        log_level=logging.DEBUG,
        storage_bucket=BUCKET,
        gem_mapper_threads=5,
    )
    # pipeline = VariantCallingPipeline(
    #     run_id="test_run",
    #     fasta_path=f"s3://{BUCKET}/fasta/hg19.fa",
    #     fasta_chunks=60,
    #     fastq_path=f"s3://{BUCKET}/fastq/SRR15068323.fastq.gz",
    #     fastq_chunks=20,
    #     log_level=logging.DEBUG,
    #     storage_bucket=BUCKET,
    #     gem_mapper_threads=5,
    # )
    # pipeline.clean_all()
    # pipeline.clean_temp_data()
    # pipeline.preprocess()
    # pipeline.alignment()
    # pipeline.reduce()
    pipeline.run_pipeline()

    # with open("hg19_alignment.json", "w") as f:
    # f.write(json.dumps(pipeline.global_stat.dump_dict()))
    with open("trityp.json", "w") as f:
        f.write(json.dumps(pipeline.global_stat.dump_dict()))
