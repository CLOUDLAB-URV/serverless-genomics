import os
import pathlib
import re
import lithops
from serverlessgenomics.parameters import PipelineParameters

class FastaPartitionerIndex:
    def __init__(self, bucket):
        self.bucket = bucket

    def __get_length(self, min_range, content, data, start_base, end_base):
        start_base -= min_range
        end_base -= min_range
        len_base = len(data[start_base:end_base].replace('\n', ''))
        # name_id offset_head offset_bases ->
        # name_id offset_head offset_bases len_bases
        content[-1] = f'{content[-1]} {len_base}'

    # Generate metadata from fasta file
    def generate_chunks(self, storage, id, key, chunk_size, obj_size, partitions):
        min_range = id * chunk_size
        max_range = int(obj_size) if id == partitions - 1 else (id + 1) * chunk_size
        data = storage.get_object(bucket=self.bucket, key=key,
                                       extra_get_args={'Range': f'bytes={min_range}-{max_range - 1}'}).decode('utf-8')
        content = []
        ini_heads = list(
            re.finditer(r"\n>", data))  # If it were '>' it would also find the ones inside the head information
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
                        # (start-1): avoid having a split sequence in the index that only has '\n'
                        match_text = list(re.finditer('.*\n', data[0:m.start()]))
                        if match_text:
                            text = match_text[0].group().split(' ')[0].replace('\n', '')
                            length_0 = len(data[match_text[0].start():m.start()].replace('\n', ''))
                            offset_0 = match_text[0].start() + min_range
                            if len(match_text) > 1:
                                offset_1 = match_text[1].start() + min_range
                                length_1 = len(data[match_text[1].start():m.start()].replace('\n', ''))
                                length_base = f"{length_0}-{length_1}"
                                offset = f"{offset_0}-{offset_1}"
                            else:
                                length_base = f"{length_0}"
                                offset = f'{offset_0}'
                            # >> offset_head offset_bases_split length/s first_line_before_space_or_\n
                            content.append(f">> <Y> {str(offset)} {length_base} ^{text}^")  # Split sequences
                        else:  # When the first header found is false, when in a split stream there is a split header that has a '>' inside (ex: >tr|...o-alpha-(1->5)-L-e...\n)
                            first_sequence = True
                            start = end = -1  # Avoid entering the following condition
                if prev != start:  # When if the current sequence base is not empty
                    if prev != -1:
                        self.__get_length(min_range, content, data, prev, start)
                    # name_id offset_head offset_bases
                    id_name = m.group().replace('\n', '').split(' ')[0].replace('>', '')
                    content.append(f"{id_name} {str(start)} {str(end)}")
                prev = end

            if len(heads) != 0 and len(ini_heads) != 0 and ini_heads[-1].start() + 1 > heads[
                -1].start():  # Check if the last head of the current one is cut. (ini_heads[-1].start() + 1): ignore '\n'
                last_seq_start = ini_heads[-1].start() + min_range + 1  # (... + 1): ignore '\n'
                self.__get_length(min_range, content, data, prev, last_seq_start)  # Add length of bases to last sequence
                text = data[last_seq_start - min_range::]
                # [<->|<_>]name_id_split offset_head
                content.append(
                    f"{'<-' if ' ' in text else '<_'}{text.split(' ')[0]} {str(last_seq_start)}")  # if '<->' there is all id
            else:  # Add length of bases to last sequence
                self.__get_length(min_range, content, data, prev, max_range)
        return content
    
    def reduce_generate_chunks(self, results):
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
                                if '<->' in list_prev[-1]:  # If the split was after a space, then there is all id
                                    name_id = param_seq_prev[0].replace('<->', '')
                                else:
                                    name_id = param_seq_prev[0].replace('<_>', '') + param[4].replace('^', '')
                                length = param[3].split('-')[1]
                                offset_head = param_seq_prev[1]
                                offset_base = param[2].split('-')[1]
                                list_prev.pop()  # Remove previous sequence
                                list_seq[0] = list_seq[0].replace(f' {param[4]}', '')  # Remove 5rt param
                                list_seq[0] = list_seq[0].replace(f' {param[2]} ',
                                                                    f' {offset_base} ')  # [offset_base_0-offset_base_1|offset_base] -> offset_base
                                list_seq[0] = list_seq[0].replace(f' {param[3]} ',
                                                                    f' {length} ')  # [length_0-length_1|length] -> length
                                list_seq[0] = list_seq[0].replace(' <Y> ', f' {offset_head} ')  # Y --> offset_head
                                list_seq[0] = list_seq[0].replace('>> ', f'{name_id} ')  # '>>' -> name_id
                            else:
                                list_seq.pop(0)
            return results

class FastaPartitionerCaller:
    def __init__(self, bucket: str):
        self.bucket = bucket
        self.storage = lithops.Storage()

    def __call__(self, my_key: str, n_workers: int, fasta_folder: str):
        # Execution   
        fexec = lithops.FunctionExecutor()

        fasta = self.storage.head_object(self.bucket, my_key)
        chunk_size = int(int(fasta['content-length']) / n_workers)
        map_iterdata = [{'key': my_key} for _ in range(n_workers)]
        funct = FastaPartitionerIndex(self.bucket)
        extra_args = {'chunk_size': chunk_size, 'obj_size': fasta['content-length'], 'partitions': n_workers}
        fexec.map_reduce(map_function=funct.generate_chunks, map_iterdata=map_iterdata, extra_args=extra_args, reduce_function=funct.reduce_generate_chunks)        
        index = fexec.get_result()        
        fexec.clean()

        self.__generate_index_file(index, fasta_folder, f'{pathlib.Path(my_key).stem}.fai') 
 
    def __generate_index_file(self, data, fasta_folder, file_name):
        data_string = ''
        for list_seq in data:
            data_string += "\n".join(list_seq) if not data_string else "\n" + "\n".join(list_seq)

        self.storage.put_object(self.bucket, f'{fasta_folder}{file_name}', str(data_string))


class FunctionsFastaIndex:
    def __init__(self, path_index_file, path_fasta_file):
        self.path_index_file = path_index_file
        self.path_fasta_file = path_fasta_file
        self.storage = lithops.Storage()  

    def get_info_sequence(self, identifier):
        length = offset_head = offset = None
        if identifier != '':
            data_index = self.storage.get_object(args.bucket, self.path_index_file).decode('utf-8').split('\n')
            size_data = len(data_index)
            i = 0
            while i < size_data:
                if identifier in data_index[i]:
                    param_seq = data_index[i].split(' ')
                    offset_head = int(param_seq[1])
                    offset = int(param_seq[2])
                    length = int(param_seq[3])
                    i+=1
                    while i < size_data and identifier in data_index[i]:
                        length += int(data_index[i].split(' ')[3])
                        i+=1
                    break
        return {'length': length, 'offset_head': offset_head, 'offset': offset}

    def get_sequences_of_range(self, min_range, max_range):
        sequences = []
        data_index = self.storage.get_object(args.bucket, self.path_index_file).decode('utf-8').split('\n')
        size_data = len(data_index)
        i = 0
        while i < size_data and int(data_index[i].split(' ')[2]) < min_range:
            i+=1

        while i < size_data and int(data_index[i].split(' ')[2]) < max_range:
            sequences.append(data_index[i].replace('\n', ''))
            i+=1
        return sequences

    def get_chunks(self, args: PipelineParameters):
        fasta_chunks = []
        fasta = self.storage.head_object(args.bucket, self.path_fasta_file)
        total_size = int(fasta['content-length'])
        fa_chunk_size = int(total_size / int(args.fasta_workers))
        data_index = self.storage.get_object(args.bucket, self.path_index_file).decode('utf-8').split('\n')
        size_data = len(data_index)
        j = 0
        min = fa_chunk_size * j
        max = fa_chunk_size * (j + 1)
        i = 0
        while max <= total_size:
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

            if i == size_data - 1 or max < int(data_index[i+1].split(' ')[1]):
                fa_chunk['last_byte+'] = max - 1 if fa_chunk_size * (j + 2) <= total_size else total_size - 1
            else:
                if max < int(data_index[i+1].split(' ')[2]):  # Split in the middle of head
                    fa_chunk['last_byte+'] = int(data_index[i+1].split(' ')[1]) - 1
                    i += 1
                elif i < size_data:
                    i += 1
                    while i + 1 < size_data and max > int(data_index[i + 1].split(' ')[1]):
                        i += 1
                    fa_chunk['last_byte+'] = max - 1
            fasta_chunks.append(fa_chunk)
            j += 1
            min = fa_chunk_size * j
            max = fa_chunk_size * (j + 1)
        return fasta_chunks

