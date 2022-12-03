from serverlessgenomics.parameters import PipelineParameters
from .fasta_partitioner_index import get_fasta_chunks
import re

def generate_alignment_iterdata(args: PipelineParameters, list_fastq: list, fasta_index: str, fasta_file_path):
    # The iterdata consists of an array where each element is a pair of a fastq chunk and a fasta chunk.
    # Since each fastq chunk needs to be paired with each fasta chunk, the total number of elements in
    # iterdata will be n_fastq_chunks * n_fasta_chunks.
    
    print("\nStarting phase: iterdata generation")
    # Get the fasta chunks
    fasta_chunks = get_fasta_chunks(fasta_index, fasta_file_path, args)

    # Generate iterdata
    num_chunks = 0
    iterdata = []
    for fastq_key in list_fastq:
        num_chunks += 1
        for i, chunk in enumerate(fasta_chunks):
            iterdata.append({'fasta_chunk': {'key_fasta': fasta_file_path, 'key_index': fasta_index, 'id': i, 'chunk': chunk}, 
                            'fastq_chunk': fastq_key, 'exec_param': args.execution_name})
    # Check iterdata
    n_fasta_chunks = len(fasta_chunks)
    if args.iterdata_n is not None:
        iterdata = iterdata[0:int(args.iterdata_n)]
        if len(iterdata) % n_fasta_chunks != 0:
            raise Exception(f"ERROR. Number of elements in iterdata must be multiple of the number of fasta chunks (max iterdata: {len(iterdata)}, data generated: {len(fasta_chunks)}).")
        else:
            num_chunks = int(re.sub('^[\s|\S]*number\':\s(\d*),[\s|\S]*$', r"\1", str(iterdata[-1]['fastq_chunk'])))
    if not iterdata:
        raise Exception('ERROR. Iterdata not generated')
    
    print("\nITERDATA LIST")
    print("   - number of fasta chunks: " + str(n_fasta_chunks))
    print("   - number of fastq chunks: " + str(len(list_fastq)))
    print("      · number chunks will be executed: " + str(num_chunks))
    print("   - fasta x fastq chunks: " + str(n_fasta_chunks * len(list_fastq)))
    print("   - number of iterdata elements: " + str(len(iterdata)))

    return iterdata, num_chunks
