import json
import matplotlib.pyplot as plt
import numpy as np
import pandas
import os

def fetch_general_data(data, dict_data):
    dict_data['pipeline_total_time'] = data['pipeline']['execution_time']
    dict_data['preprocessing_total_time'] = data['pipeline']['preprocess_phase']['preprocess']['execution_time']
    dict_data['map_total_time'] = data['pipeline']['alignReads_phase']['align_reads']['execution_time']
    dict_data['reduce_total_time'] = data['pipeline']['reduce_phase']['reduce']['execution_time']
    return dict_data


def plot_map_one(data: dict):
    function_details = data['pipeline']['alignReads_phase']['align_reads']['phases']['gem_indexer_mapper']['function_details']
    
    i = 0
    download_fastq = []
    download_fasta = []
    gem_indexer = []
    map_index_and_filter_map = []
    compress_index = []
    compress_map = []
    upload_index = []
    upload_map = []
    end = []
    
    for elem in function_details:
        keys = elem.keys()
        for k in keys:
            timestamps = elem[k]['timestamps']
            start = timestamps['start']
            download_fastq.append(timestamps['download_fastq'] - start)
            download_fasta.append(timestamps['download_fasta'] - timestamps['download_fastq'])
            gem_indexer.append(timestamps['gem_indexer'] - timestamps['download_fasta'])
            map_index_and_filter_map.append(timestamps['map_index_and_filter_map'] - timestamps['gem_indexer'])
            compress_index.append(timestamps['compress_index'] - timestamps['map_index_and_filter_map'])
            compress_map.append(timestamps['compress_map'] - timestamps['compress_index'])
            upload_index.append(timestamps['upload_index'] - timestamps['compress_map'])
            upload_map.append(timestamps['upload_map'] - timestamps['upload_index'])
            end.append(timestamps['end'] - timestamps['upload_map'])
        i += 1
            
    df = pandas.DataFrame({
        'download_fastq': download_fasta,
        'download_fasta': gem_indexer,
        'gem_indexer': map_index_and_filter_map,
        'map_index_and_filter_map': compress_index,
        'compress_index': compress_map,
        'compress_map': upload_index,
        'upload_index': upload_map,
        'upload_map': end
    })
    
    ax = df.plot.barh(stacked=True, title='Map Phase One')
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
    
    ax.set_yticks([0, i])
    ax.set_yticklabels([0, i])
    
    ax.set_xlabel("time in seconds")
    ax.set_ylabel("functions")
    
    fig = ax.get_figure()
    fig.savefig("./stats/map_phase_one.png", bbox_inches='tight', dpi=400)


def plot_map_two(data: dict):
    function_details = data['pipeline']['alignReads_phase']['align_reads']['phases']['filter_index_to_mpileup']['function_details']
    
    i = 0
    download_fasta_chunk = []
    download_map_file = []
    download_index = []
    map_file_index_correction = []
    gempileup_run = []
    upload_mpileup = []
    end = []
    
    for elem in function_details:
        keys = elem.keys()
        for k in keys:
            timestamps = elem[k]['timestamps']
            start = timestamps['start']
            download_fasta_chunk.append(timestamps['download_fasta_chunk'] - start)
            download_map_file.append(timestamps['download_map_file'] - timestamps['download_fasta_chunk'])
            download_index.append(timestamps['download_index'] - timestamps['download_map_file'])
            map_file_index_correction.append(timestamps['map_file_index_correction'] - timestamps['download_index'])
            gempileup_run.append(timestamps['gempileup_run'] - timestamps['map_file_index_correction'])
            upload_mpileup.append(timestamps['upload_mpileup'] - timestamps['gempileup_run'])
            end.append(timestamps['end'] - timestamps['upload_mpileup'])
        i += 1
            
    df = pandas.DataFrame({
        'download_fasta_chunk': download_map_file,
        'download_map_file': download_map_file,
        'download_map_file': download_index,
        'download_index': map_file_index_correction,
        'map_file_index_correction': gempileup_run,
        'gempileup_run': upload_mpileup,
        'upload_mpileup': end
    })
    
    ax = df.plot.barh(stacked=True, title='Map Phase Two')
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
    
    ax.set_yticks([0, i])
    ax.set_yticklabels([0, i])
    
    ax.set_xlabel("time in seconds")
    ax.set_ylabel("functions")
    
    fig = ax.get_figure()
    fig.savefig("./stats/map_phase_two.png", bbox_inches='tight', dpi=400)


def plot_index_correction(data: dict):
    function_details = data['pipeline']['alignReads_phase']['align_reads']['phases']['index_correction']['function_details']
    
    i = 0
    download_indexes = []
    merge_gem = []
    filter_merged = []
    compress_corrected_index = []
    upload_corrected_index = []
    end = []
    
    for elem in function_details:
        keys = elem.keys()
        for k in keys:
            timestamps = elem[k]['timestamps']
            start = timestamps['start']
            download_indexes.append(timestamps['download_indexes'] - start)
            merge_gem.append(timestamps['merge_gem'] - timestamps['download_indexes'])
            filter_merged.append(timestamps['filter_merged'] - timestamps['merge_gem'])
            compress_corrected_index.append(timestamps['compress_corrected_index'] - timestamps['filter_merged'])
            upload_corrected_index.append(timestamps['upload_corrected_index'] - timestamps['compress_corrected_index'])
            end.append(timestamps['end'] - timestamps['upload_corrected_index'])
        i += 1
            
    df = pandas.DataFrame({
        'download_indexes': merge_gem,
        'merge_gem': filter_merged,
        'filter_merged': compress_corrected_index,
        'compress_corrected_index': upload_corrected_index,
        'upload_corrected_index': end
    })
    
    ax = df.plot.barh(stacked=True, title='Index Correction')
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
    
    ax.set_yticks([0, i])
    ax.set_yticklabels([0, i])
    
    ax.set_xlabel("time in seconds")
    ax.set_ylabel("functions")
    
    fig = ax.get_figure()
    fig.savefig("./stats/index_correction.png", bbox_inches='tight', dpi=400)
    

def plot_distribute_indexes(data: dict):
    function_details = data['pipeline']['reduce_phase']['reduce']['phases']['distribute_indexes']['function_details']
    
    i = 0
    s3_queries = []
    distribute_indexes = []
    end = []
    
    for elem in function_details:
        k = list(elem.keys())[0]
        timestamps = elem[k]['timestamps']
        start = timestamps['start']
        s3_queries.append(timestamps['s3_queries'] - start)
        distribute_indexes.append(timestamps['distribute_indexes'] - timestamps['s3_queries'])
        end.append(timestamps['end'] - timestamps['distribute_indexes'])
        i += 1
            
    df = pandas.DataFrame({
        's3_queries': distribute_indexes,
        'distribute_indexes': end
    })
    
    ax = df.plot.barh(stacked=True, title='Distribute Indexes')
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
    
    ax.set_yticks([0, i])
    ax.set_yticklabels([0, i])
    
    ax.set_xlabel("time in seconds")
    ax.set_ylabel("functions")
    
    fig = ax.get_figure()
    fig.savefig("./stats/distribute_indexes.png", bbox_inches='tight', dpi=400)


def plot_reduce(data: dict):
    function_details = data['pipeline']['reduce_phase']['reduce']['phases']['reduce_function']['function_details']
    
    i = 0
    s3_queries = []
    mpileup_merge_reduce = []
    upload_part = []
    end = []
    
    for elem in function_details:
        keys = elem.keys()
        for k in keys:
            timestamps = elem[k]['timestamps']
            start = timestamps['start']
            s3_queries.append(timestamps['s3_queries'] - start)
            mpileup_merge_reduce.append(timestamps['mpileup_merge_reduce'] - timestamps['s3_queries'])
            upload_part.append(timestamps['upload_part'] - timestamps['mpileup_merge_reduce'])
            end.append(timestamps['end'] - timestamps['upload_part'])
        i += 1
            
    df = pandas.DataFrame({
        's3_queries': mpileup_merge_reduce,
        'mpileup_merge_reduce': upload_part,
        'upload_part': end
    })
    
    ax = df.plot.barh(stacked=True, title='Reduce Function')
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
    
    ax.set_yticks([0, i])
    ax.set_yticklabels([0, i])
    
    ax.set_xlabel("time in seconds")
    ax.set_ylabel("functions")
    
    fig = ax.get_figure()
    fig.savefig("./stats/reduce_functions.png", bbox_inches='tight', dpi=400)


if __name__ == '__main__':
    with open("/home/agabriel/Downloads/logs_stats.json") as json_read:
        data: dict = json.load(json_read)
    
    if not os.path.exists('stats'):
        os.makedirs('stats')
    
    plot_map_one(data)
    plot_map_two(data)
    plot_index_correction(data)
    plot_distribute_indexes(data)
    plot_reduce(data)