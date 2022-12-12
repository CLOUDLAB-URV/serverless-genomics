import os
import pathlib
import re
import lithops
from serverlessgenomics.parameters import PipelineParameters

class FastaPartitionerIndex:
    def __init__(self, bucket):
        self.bucket = bucket

    # Generate metadata from fasta file
    def generate_chunks(self, storage, id, key, chunk_size, obj_size, partitions):
        '''
        Searches for the different sequences in the given chunk and returns a list of them
        '''
        min_range = id * chunk_size
        max_range = int(obj_size) if id == partitions - 1 else (id + 1) * chunk_size
        data = storage.get_object(bucket=self.bucket, key=key, extra_get_args={'Range': f'bytes={min_range}-{max_range - 1}'}).decode('utf-8')
        content = []
        ini_heads = list(re.finditer(r"\n>", data))  # If it were '>' it would also find the ones inside the head information
        heads = list(re.finditer(r">.+\n", data))

        if ini_heads or data[0] == '>':  # If the list is not empty or there is > in the first byte
            first_sequence = True
            prev = -1
            for m in heads:
                start = min_range + m.start()
                end = min_range + m.end()
                if first_sequence:
                    first_sequence = False
                    if id > 0 and start - 1 > min_range:  # If it is not the worker of the first part of the file and in addition it
                        # turns out that the partition begins in the middle of the base of a sequence.
                        # (start-1): avoid having a split sequence in the index that only has '\n'.
                        match_text = list(re.finditer('.*\n', data[0:m.start()]))
                        if match_text and len(match_text) > 1:
                            text = match_text[0].group().split(' ')[0].replace('\n', '')
                            offset = match_text[1].start() + min_range
                            # >> offset_head offset_bases_split ^first_line_before_space_or_\n^
                            content.append(f">> <Y> {str(offset)} ^{text}^")  # Split sequences
                        else:  # When the first header found is false, when in a split stream there is a split header that has a '>' inside (ex: >tr|...o-alpha-(1->5)-L-e...\n)
                            first_sequence = True
                if prev != start:  # When if the current sequence base is not empty                    
                    id_name = m.group().replace('\n', '').split(' ')[0].replace('>', '')
                    # name_id offset_head offset_bases
                    content.append(f"{id_name} {str(start)} {str(end)}")
                prev = end

            if len(heads) != 0 and len(ini_heads) != 0 and ini_heads[-1].start() + 1 > heads[
                -1].start():  # Check if the last head of the current one is cut. (ini_heads[-1].start() + 1): ignore '\n'
                last_seq_start = ini_heads[-1].start() + min_range + 1  # (... + 1): ignore '\n'
                text = data[last_seq_start - min_range::]
                # [<->|<_>]name_id_split offset_head
                content.append(f"{'<-' if ' ' in text else '<_'}{text.split(' ')[0]} {str(last_seq_start)}")  # if '<->' there is all id
        return content
    
    def reduce_generate_chunks(self, results):
        '''
        Check the first and last position of each list and solve the sequences that were split. Return a list of lists
        '''
        if len(results) > 1:            
            results = list(filter(None, results))
            for i, list_seq in enumerate(results):
                if i > 0:
                    list_prev = results[i - 1]
                    if list_prev and list_seq:  # If it is not empty the current and previous dictionary
                        param = list_seq[0].split(' ')
                        seq_prev = list_prev[-1]
                        param_seq_prev = seq_prev.split(' ')
                        if '>>' in list_seq[0]:  # If the first sequence is split
                            if '<->' in seq_prev or '<_>' in seq_prev:
                                if '<->' in seq_prev:  # If the split was after a space, then there is all id
                                    name_id = param_seq_prev[0].replace('<->', '')
                                else:
                                    name_id = param_seq_prev[0].replace('<_>', '') + param[3].replace('^', '')                       
                                list_seq[0] = self.__rename_seq(list_seq[0], param, name_id, param_seq_prev[1], param[2])
                            else:
                                list_seq[0] = seq_prev
                            list_prev.pop()  # Remove previous sequence
        return results

    def __rename_seq(sequence, param, name_id, offset_head, offset_base):    
        sequence = sequence.replace(f' {param[3]}', '')  # Remove 3rt param
        sequence = sequence.replace(f' {param[2]} ', f' {offset_base} ')  # offset_base -> offset_base
        sequence = sequence.replace(' <Y> ', f' {offset_head} ')  # Y --> offset_head
        sequence = sequence.replace('>> ', f'{name_id} ')  # '>>' -> name_id
        return sequence

    def send_index_file_to_storage(self, data, fasta_folder, file_name, storage):
        '''
        Function that loads the generated index file into the storage. Read a list of lists and converts it into a string, 
        then uploads it to the storage
        '''
        data_string = ''
        for list_seq in data:
            data_string += "\n".join(list_seq) if not data_string else "\n" + "\n".join(list_seq)

        storage.put_object(self.bucket, f'{fasta_folder}{file_name}', str(data_string))


def fasta_partitioner_caller(bucket: str, my_key: str, n_workers: int, fasta_folder: str): 
    '''
    Call function of the index generation functions
    '''
    fexec = lithops.FunctionExecutor()
    storage = lithops.Storage()

    # Get metaddata of the fasta file
    fasta_length = storage.head_object(bucket, my_key)['content-length']
    chunk_size = int(int(fasta_length) / n_workers)
    funct = FastaPartitionerIndex(bucket)
  
    # Preparing parameters for serverless functions 
    map_iterdata = [{'key': my_key}] * n_workers
    extra_args = {'chunk_size': chunk_size, 'obj_size': fasta_length, 'partitions': n_workers}
    
    fexec.map_reduce(map_function=funct.generate_chunks, map_iterdata=map_iterdata, extra_args=extra_args, reduce_function=funct.reduce_generate_chunks)        
    index = fexec.get_result()        
    fexec.clean()

    if not index:  # Check if the data was generated 
        raise Exception('ERROR: the generated index file is empty.')

    funct.send_index_file_to_storage(index, fasta_folder, f'{pathlib.Path(my_key).stem}.fai', storage) 
    
    

def get_fasta_chunks(path_index_file: str, path_fasta_file: str, args: PipelineParameters):
    '''
    Read the index file and get the chunks based on the number of workers set
    '''
    storage = lithops.Storage() 

    fasta_chunks = []
    total_size = int(storage.head_object(args.bucket, path_fasta_file)['content-length'])
    fa_chunk_size = int(total_size / int(args.fasta_workers))
    data_index = storage.get_object(args.bucket, path_index_file).decode('utf-8').split('\n')      
    size_data = len(data_index)
    
    i = j = 0
    min = fa_chunk_size * j
    max = fa_chunk_size * (j + 1)
    while max <= total_size:
        # Find first full/half sequence of the chunk
        if int(data_index[i].split(' ')[1]) <= min < int(data_index[i].split(' ')[2]):   # In the head
            fa_chunk = {'offset_head': int(data_index[i].split(' ')[1]), 'offset_base': int(data_index[i].split(' ')[2])}
        elif i == size_data - 1 or min < int(data_index[i + 1].split(' ')[1]):  # In the base
            fa_chunk = {'offset_head': int(data_index[i].split(' ')[1]), 'offset_base': min}
        elif i < size_data:
            i += 1
            while i + 1 < size_data and min > int(data_index[i + 1].split(' ')[1]):
                i += 1
            if min < int(data_index[i].split(' ')[2]):
                fa_chunk = {'offset_head': int(data_index[i].split(' ')[1]), 'offset_base': int(data_index[i].split(' ')[2])}
            else:
                fa_chunk = {'offset_head': int(data_index[i].split(' ')[1]), 'offset_base': min}
        else:
            raise Exception('ERROR: there was a problem getting the first byte of a fasta chunk.')
        # Find last full/half sequence of the chunk
        if i == size_data - 1 or max < int(data_index[i + 1].split(' ')[1]):
            fa_chunk['last_byte+'] = max - 1 if fa_chunk_size * (j + 2) <= total_size else total_size - 1
        else:
            if max < int(data_index[i + 1].split(' ')[2]):  # Split in the middle of head
                fa_chunk['last_byte+'] = int(data_index[i + 1].split(' ')[1]) - 1
                i += 1
            elif i < size_data:
                i += 1
                while i + 1 < size_data and max > int(data_index[i + 1].split(' ')[1]):
                    i += 1
                fa_chunk['last_byte+'] = max - 1
            else:
                raise Exception('ERROR: there was a problem getting the last byte of a fasta chunk.')
                
        fasta_chunks.append(fa_chunk)
        j += 1
        min = fa_chunk_size * j
        max = fa_chunk_size * (j + 1)
    return fasta_chunks


