import os
import subprocess
import xml
import logging

import requests

from serverlessgenomics.pipelineparams import PipelineParameters

logger = logging.getLogger(__name__)

SRA_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def get_sra_metadata(pipeline_params: PipelineParameters) -> int:
    params = {"db": "sra", "id": pipeline_params.sra_accession, "retmode": "xml"}

    response = requests.get(SRA_EFETCH_URL, params=params)

    if response.status_code == 200:
        xml_data = response.text
        root = xml.etree.ElementTree.fromstring(xml_data)

        for run in root.iter("RUN"):
            reads = int(run.get("total_spots"))
            logger.debug("Read %d total reads from efetch for sequence %s", reads, pipeline_params.sra_accession)
            return reads
    else:
        raise Exception(f"Error fetching metadata for {pipeline_params.sra_accession}: {response.status_code}")


def fetch_fastq_chunk_sra(seq_name: str, fastq_chunk: dict, target_filename: str):
    """
    Function to retrieve the relevant SRA chunk using fastq-dump and save it to object storage
    """

    start_read = int(fastq_chunk["read_0"])
    end_read = int(fastq_chunk["read_1"])

    # To suppress a warning that appears the first time vdb-config is used
    proc = subprocess.run(["vdb-config", "-i"], check=True, capture_output=True, text=True)
    print(proc.stdout)
    print(proc.stderr)
    # Report cloud identity so it can take data from s3 needed to be executed only once per vm
    proc = subprocess.run(["vdb-config", "--report-cloud-identity", "yes"], check=True, capture_output=True, text=True)
    print(proc.stdout)
    print(proc.stderr)

    # Run fastq-dump with the specified range of reads, splits files in two files if paired end
    proc = subprocess.run(
        ["fastq-dump", "--split-files", seq_name, "-X", str(start_read), "-N", str(end_read)], check=True,
        capture_output=True, text=True
    )
    print(proc.stdout)
    print(proc.stderr)

    fastqdump_output_filename = f"{seq_name}_1.fastq"

    os.rename(fastqdump_output_filename, target_filename)

    # Save the output file to the desired target filename in object storage

    print(f"Finished fetching chunk {fastq_chunk['chunk_id']} and saved to {target_filename}")
