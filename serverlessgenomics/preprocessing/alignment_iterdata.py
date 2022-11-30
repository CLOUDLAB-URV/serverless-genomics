from serverlessgenomics.parameters import PipelineParameters
from .fasta_partitioner_index import FunctionsFastaIndex
import re

class GenerateAlignmentIterdata:
    def __init__(self, args: PipelineParameters, list_fastq: list, fasta_index: str, fasta_file_path):
        self.args = args
        self.list_fastq = list_fastq
        self.fasta_index = fasta_index
        self.fasta_file_path = fasta_file_path

    def __call__(self):
        """
        Creates the lithops iterdata from the fasta and fastq chunk lists
        """
        # The iterdata consists of an array where each element is a pair of a fastq chunk and a fasta chunk.
        # Since each fastq chunk needs to be paired with each fasta chunk, the total number of elements in
        # iterdata will be n_fastq_chunks * n_fasta_chunks.
        
        print("\nStarting phase: iterdata generation")
        functions = FunctionsFastaIndex(self.fasta_index, self.fasta_file_path)
        fasta_chunks = functions.get_chunks(self.args)

        iterdata, num_chunks = self.__generate_iterdata(fasta_chunks)

        return self.__check_iterdata(iterdata, fasta_chunks, num_chunks)

    def __generate_iterdata(self, fasta_chunks):
        num_chunks = 0
        iterdata = []
        for fastq_key in self.list_fastq:
            num_chunks += 1
            for i, chunk in enumerate(fasta_chunks):
                iterdata.append({'fasta_chunk': {'key_fasta': self.fasta_file_path, 'key_index': self.fasta_index, 'id': i, 'chunk': chunk}, 
                                'fastq_chunk': fastq_key, 'exec_param': self.args.execution_name})

        return iterdata, num_chunks

    def __check_iterdata(self, iterdata, fasta_chunks, num_chunks):
        n_fasta_chunks = len(fasta_chunks)
        n_fastq = len(self.list_fastq)
        
        if self.args.iterdata_n is not None:
            iterdata = iterdata[0:int(self.args.iterdata_n)]
            if len(iterdata) % n_fasta_chunks != 0:
                raise Exception(
                    f"ERROR. Number of elements in iterdata must be multiple of the number of fasta chunks (max iterdata: {len(iterdata)}, data generated: {len(fasta_chunks)}).")
            else:
                num_chunks = int(re.sub('^[\s|\S]*number\':\s(\d*),[\s|\S]*$', r"\1", str(iterdata[-1]['fastq_chunk'])))
        if not iterdata:
            raise Exception('ERROR. Iterdata not generated')

        print("\nITERDATA LIST")
        print("   - number of fasta chunks: " + str(n_fasta_chunks))
        print("   - number of fastq chunks: " + str(n_fastq))
        print("      · number chunks will be executed: " + str(num_chunks))
        print("   - fasta x fastq chunks: " + str(n_fasta_chunks * n_fastq))
        print("   - number of iterdata elements: " + str(len(iterdata)))

        return iterdata, num_chunks