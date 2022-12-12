from __future__ import annotations

import math
import os
import re
import logging
import itertools
from typing import TYPE_CHECKING
from functools import reduce

from ..utils import try_head_object
import bz2

if TYPE_CHECKING:
    from serverlessgenomics.parameters import PipelineRun, Lithops

logger = logging.getLogger(__name__)


def create_index_chunked(storage, id, fasta_path, chunk_size, fasta_size, num_chunks):
    """
    Lithops callee function (map)
    Generate partial index of a chunk of a FASTA file
    """
    min_range = id * chunk_size
    max_range = int(fasta_size) if id == num_chunks - 1 else (id + 1) * chunk_size
    data = storage.get_object(bucket=fasta_path.bucket, key=fasta_path.key,
                              extra_get_args={'Range': f'bytes={min_range}-{max_range - 1}'}).decode('utf-8')
    content = []
    # If it were '>' it would also find the ones inside the head information
    ini_heads = list(re.finditer(r"\n>", data))
    heads = list(re.finditer(r">.+\n", data))

    if ini_heads or data[0] == '>':  # If the list is not empty or there is > in the first byte
        first_sequence = True
        prev = -1
        for m in heads:
            start = min_range + m.start()
            end = min_range + m.end()
            if first_sequence:
                first_sequence = False
                if id > 0 and start - 1 > min_range:
                    # If it is not the worker of the first part of the file and in addition it
                    # turns out that the partition begins in the middle of the base of a sequence.
                    # (start-1): avoid having a split sequence in the index that only has '\n'.
                    match_text = list(re.finditer('.*\n', data[0:m.start()]))
                    if match_text and len(match_text) > 1:
                        text = match_text[0].group().split(' ')[0].replace('\n', '')
                        offset = match_text[1].start() + min_range
                        # >> offset_head offset_bases_split ^first_line_before_space_or_\n^
                        content.append(f">> <Y> {str(offset)} ^{text}^")  # Split sequences
                    else:
                        # When the first header found is false, when in a split stream there is a split header
                        # that has a '>' inside (ex: >tr|...o-alpha-(1->5)-L-e...\n)
                        first_sequence = True
            if prev != start:  # When if the current sequence base is not empty
                # name_id offset_head offset_bases
                id_name = m.group().replace('\n', '').split(' ')[0].replace('>', '')
                content.append(f"{id_name} {str(start)} {str(end)}")
            prev = end

        # Check if the last head of the current one is cut. (ini_heads[-1].start() + 1): ignore '\n'
        if len(heads) != 0 and len(ini_heads) != 0 and ini_heads[-1].start() + 1 > heads[-1].start():
            last_seq_start = ini_heads[-1].start() + min_range + 1  # (... + 1): ignore '\n'
            text = data[last_seq_start - min_range::]
            # [<->|<_>]name_id_split offset_head
            # if '<->' there is all id
            content.append(f"{'<-' if ' ' in text else '<_'}{text.split(' ')[0]} {str(last_seq_start)}")

    return content


def reduce_chunked_indexes(results, storage):
    """
    Lithops callee function (reduce)
    Reduce partial indexes of chunked FASTA file into a single index for the entire file
    """
    # TODO replace with reduce extra args (pending issue lithops PR)
    bucket = os.getenv('BUCKET')
    faidx_key = os.getenv('FAIDX_KEY')

    if len(results) > 1:
        results = list(filter(None, results))
        for i, list_seq in enumerate(results):
            if i > 0:
                list_prev = results[i - 1]
                # If it is not empty the current and previous dictionary
                if list_prev and list_seq:
                    param = list_seq[0].split(' ')
                    seq_prev = list_prev[-1]
                    param_seq_prev = seq_prev.split(' ')

                    # If the first sequence is split
                    if '>>' in list_seq[0]:
                        if '<->' in seq_prev or '<_>' in seq_prev:
                            # If the split was after a space, then there is all id
                            if '<->' in seq_prev:
                                name_id = param_seq_prev[0].replace('<->', '')
                            else:
                                name_id = param_seq_prev[0].replace('<_>', '') + param[3].replace('^', '')
                            list_seq[0] = rename_sequence(list_seq[0], param, name_id, param_seq_prev[1], param[2])
                        else:
                            list_seq[0] = seq_prev
                        # Remove previous sequence
                        list_prev.pop()

    num_sequences = reduce(lambda x, y: x + y, map(lambda r: len(r), results))

    index_result = bz2.compress(b'\n'.join(s.encode('utf-8') for s in itertools.chain(*results)))
    s3 = storage.get_client()
    s3.put_object(Bucket=bucket, Key=faidx_key, Body=index_result, Metadata={'num_sequences': str(num_sequences)})
    return num_sequences


def rename_sequence(sequence, param, name_id, offset_head, offset_base):
    sequence = sequence.replace(f' {param[3]}', '')  # Remove 3rt param
    sequence = sequence.replace(f' {param[2]} ', f' {offset_base} ')  # offset_base -> offset_base
    sequence = sequence.replace(' <Y> ', f' {offset_head} ')  # Y --> offset_head
    sequence = sequence.replace('>> ', f'{name_id} ')  # '>>' -> name_id
    return sequence


def generate_faidx_from_s3(pipeline_params: PipelineRun, lithops: Lithops):
    fasta_head = try_head_object(lithops.storage, pipeline_params.fasta_path.bucket, pipeline_params.fasta_path.key)
    if fasta_head is None:
        raise Exception(f'fasta file with key {pipeline_params.fastq_path} does not exists')

    faidx_head = try_head_object(lithops.storage, pipeline_params.storage_bucket, pipeline_params.faidx_key)
    if faidx_head is not None:
        logger.debug('Faidx for %s found', pipeline_params.fasta_path.stem)
        num_sequences = int(faidx_head['x-amz-meta-num_sequences'])
    else:
        logger.info('Faidx for %s not found, generating fasta index file', pipeline_params.fasta_path.stem)

        fasta_head = lithops.storage.head_object(pipeline_params.fasta_path.bucket, pipeline_params.fasta_path.key)
        fasta_file_sz = int(fasta_head['content-length'])
        chunk_size = math.ceil(fasta_file_sz / pipeline_params.fasta_chunks)

        map_iterdata = [{'fasta_path': pipeline_params.fasta_path} for _ in range(pipeline_params.fasta_chunks)]
        extra_args = {'chunk_size': chunk_size,
                      'fasta_size': fasta_file_sz,
                      'num_chunks': pipeline_params.fasta_chunks}
        extra_env = {'BUCKET': pipeline_params.storage_bucket, 'FAIDX_KEY': pipeline_params.faidx_key}
        num_sequences = lithops.invoker.map_reduce(map_function=create_index_chunked,
                                                   map_iterdata=map_iterdata,
                                                   extra_args=extra_args, extra_env=extra_env,
                                                   reduce_function=reduce_chunked_indexes)

        logger.info('Generated faidx for FASTA %s (read %d sequences)', pipeline_params.fasta_path.stem, num_sequences)

    logger.info('Read %d sequences from FASTA %s', num_sequences, pipeline_params.fasta_path.stem)
    return num_sequences


def get_fasta_byte_ranges(pipeline_params: PipelineRun, lithops: Lithops, num_sequences):
    fasta_chunks = []
    fasta_file_head = lithops.storage.head_object(pipeline_params.fasta_path.bucket, pipeline_params.fasta_path.key)
    fasta_file_sz = int(fasta_file_head['content-length'])
    fa_chunk_size = int(fasta_file_sz / int(pipeline_params.fasta_chunks))
    compressed_faidx = lithops.storage.get_object(pipeline_params.storage_bucket, pipeline_params.faidx_key)
    faidx = bz2.decompress(compressed_faidx).decode('utf-8').split('\n')

    i = j = 0
    min = fa_chunk_size * j
    max = fa_chunk_size * (j + 1)
    while max <= fasta_file_sz:
        # Find first full/half sequence of the chunk
        if int(faidx[i].split(' ')[1]) <= min < int(faidx[i].split(' ')[2]):  # In the head
            fa_chunk = {'offset_head': int(faidx[i].split(' ')[1]),
                        'offset_base': int(faidx[i].split(' ')[2])}
        elif i == num_sequences - 1 or min < int(faidx[i + 1].split(' ')[1]):  # In the base
            fa_chunk = {'offset_head': int(faidx[i].split(' ')[1]), 'offset_base': min}
        elif i < num_sequences:
            i += 1
            while i + 1 < num_sequences and min > int(faidx[i + 1].split(' ')[1]):
                i += 1
            if min < int(faidx[i].split(' ')[2]):
                fa_chunk = {'offset_head': int(faidx[i].split(' ')[1]),
                            'offset_base': int(faidx[i].split(' ')[2])}
            else:
                fa_chunk = {'offset_head': int(faidx[i].split(' ')[1]), 'offset_base': min}
        else:
            raise Exception('ERROR: there was a problem getting the first byte of a fasta chunk.')
        # Find last full/half sequence of the chunk
        if i == num_sequences - 1 or max < int(faidx[i + 1].split(' ')[1]):
            fa_chunk['last_byte'] = max - 1 if fa_chunk_size * (j + 2) <= fasta_file_sz else fasta_file_sz - 1
        else:
            if max < int(faidx[i + 1].split(' ')[2]):  # Split in the middle of head
                fa_chunk['last_byte'] = int(faidx[i + 1].split(' ')[1]) - 1
                i += 1
            elif i < num_sequences:
                i += 1
                while i + 1 < num_sequences and max > int(faidx[i + 1].split(' ')[1]):
                    i += 1
                fa_chunk['last_byte'] = max - 1
            else:
                raise Exception('ERROR: there was a problem getting the last byte of a fasta chunk.')

        fa_chunk['chunk_id'] = j
        fasta_chunks.append(fa_chunk)
        j += 1
        min = fa_chunk_size * j
        max = fa_chunk_size * (j + 1)

    return fasta_chunks


def prepare_fasta_chunks(pipeline_params: PipelineRun, lithops: Lithops):
    """
    Calculate fasta byte ranges and metadata for chunks of a pipeline run, generate faidx index if needed
    """
    # Get number of sequences from fasta file, generate faidx file if needed
    num_sequences = generate_faidx_from_s3(pipeline_params, lithops)
    fasta_chunks = get_fasta_byte_ranges(pipeline_params, lithops, num_sequences)
    return fasta_chunks
