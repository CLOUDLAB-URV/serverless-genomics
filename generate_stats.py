import json
import matplotlib.pyplot as plt
import numpy as np

def fetch_general_data(data, dict_data):
    #General data
    dict_data['pipeline_total_time'] = data['pipeline']['execution_time']
    dict_data['preprocessing_total_time'] = data['pipeline']['preprocess_phase']['preprocess']['execution_time']
    dict_data['map_total_time'] = data['pipeline']['alignReads_phase']['align_reads']['execution_time']
    dict_data['reduce_total_time'] = data['pipeline']['reduce_phase']['reduce']['execution_time']
    return dict_data

def fetch_map_data(data, dict_data):
    #Map: Stage 1 Info
    exec_times = []
    data_transfers = []
    for result in data['pipeline']['alignReads_phase']['align_reads']['subprocesses']['gem_indexer_mapper']['subprocesses']:
        keys = list(result.keys())
        exec_times.append(result[keys[0]]['execution_time'])
        data_transfers.append(result[keys[1]])
    dict_data['map1_time'] = np.average(exec_times)
    
    fastq_fetch = []
    fasta_fetch = []
    upload_index = []
    upload_map = []
    for datat in data_transfers:
        fastq_fetch.append(datat['fastq_fetch']['execution_time'])
        fasta_fetch.append(datat['fasta_fetch']['execution_time'])
        upload_index.append(datat['upload_index']['execution_time'])
        upload_map.append(datat['upload_map']['execution_time'])
    dict_data["map1_fastq_fetch"] = np.average(fastq_fetch)
    dict_data["map1_fasta_fetch"] = np.average(fasta_fetch)
    dict_data["map1_upload_index"] = np.average(upload_index)
    dict_data["map1_upload_map"] = np.average(upload_map)
    
    
    #Map: Index Correction Info
    exec_times = []
    data_transfers = []
    for result in data['pipeline']['alignReads_phase']['align_reads']['subprocesses']['index_correction']['subprocesses']:
        keys = list(result.keys())
        exec_times.append(result[keys[0]]['execution_time'])
        data_transfers.append(result[keys[1]])
    dict_data['index_correction_time'] = np.average(exec_times)
    
    download_indexes = []
    upload_corrected_index = []
    for datat in data_transfers:
        download_indexes.append(datat['download_indexes']['execution_time'])
        upload_corrected_index.append(datat['upload_corrected_index']['execution_time'])
    dict_data["download_indexes"] = np.average(download_indexes)
    dict_data["upload_corrected_index"] = np.average(upload_corrected_index)
    
    
    #Map: Stage 2 Info
    exec_times = []
    data_transfers = []
    for result in data['pipeline']['alignReads_phase']['align_reads']['subprocesses']['filter_index_to_mpileup']['subprocesses']:
        keys = list(result.keys())
        exec_times.append(result[keys[0]]['execution_time'])
        data_transfers.append(result[keys[1]])
    dict_data['map2_time'] = np.average(exec_times)
    
    download_fasta_chunk = []
    download_map = []
    download_corrected_index = []
    upload_mpileup = []
    for datat in data_transfers:
        download_fasta_chunk.append(datat['download_fasta_chunk']['execution_time'])
        download_map.append(datat['download_map']['execution_time'])
        download_corrected_index.append(datat['download_corrected_index']['execution_time'])
        upload_mpileup.append(datat['upload_mpileup']['execution_time'])
    dict_data["map2_download_fasta_chunk"] = np.average(download_fasta_chunk)
    dict_data["map2_download_map"] = np.average(download_map)
    dict_data["map2_download_corrected_index"] = np.average(download_corrected_index)
    dict_data["map2_upload_mpileup"] = np.average(upload_mpileup)

    return dict_data


def fetch_reduce_data(data, dict_data):
    #Distribute indexes
    exec_times = []
    data_transfers = []
    for result in data['pipeline']['reduce_phase']['reduce']['subprocesses']['distribute_indexes']['subprocesses']:
        keys = list(result.keys())
        exec_times.append(result[keys[0]]['execution_time'])
        data_transfers.append(result[keys[1]]['s3_select']['execution_time'])
    dict_data['distribute_indexes_time'] = np.average(exec_times)
    dict_data["distribute_indexes_select"] = np.average(data_transfers)
    
    
    #Reduce function
    exec_times = []
    data_transfers = []
    for result in data['pipeline']['reduce_phase']['reduce']['subprocesses']['reduce_function']['subprocesses']:
        keys = list(result.keys())
        exec_times.append(result[keys[0]]['execution_time'])
        data_transfers.append(result[keys[1]])
    dict_data['reduce_function_time'] = np.average(exec_times)
    
    s3_select = []
    upload_part = []
    for datat in data_transfers:
        s3_select.append(datat['s3_select']['execution_time'])
        upload_part.append(datat['upload_part']['execution_time'])
    dict_data["reduce_func_select"] = np.average(s3_select)
    dict_data["reduce_func_upload"] = np.average(upload_part)
    
    
    #Final merge
    download_part = []
    upload_part = []
    for result in data['pipeline']['reduce_phase']['reduce']['subprocesses']['final_merge']['subprocesses']:
        keys = list(result.keys())
        download_part.append(result[keys[0]]['execution_time'])
        upload_part.append(result[keys[1]]['execution_time'])
    dict_data['final_merge_download'] = np.average(download_part)
    dict_data['final_merge_upload'] = np.average(upload_part)
    dict_data['final_merge_total'] = dict_data['final_merge_download'] + dict_data['final_merge_upload']
    
    return dict_data


def f2(num):
    return f"{num:.2f}"


def display_general(dict_data):
    print("***GENERAL INFORMATION***")
    print("Pipeline execution time: " + f2(dict_data['pipeline_total_time']) + "s")
    print("Map execution time: " + f2(dict_data['preprocessing_total_time']) + "s")
    print("Reduce execution time: " + f2(dict_data['reduce_total_time']) + "s")
    print("")
    

def display_map(dict_data):
    print("**MAP PHASE**")
    print("*MAP STAGE 1")
    print("Average time: " + f2(dict_data['map1_time']) + "s")
    total_time = float(dict_data['map1_fastq_fetch']) + float(dict_data['map1_fasta_fetch'])
    perc = f2(total_time / dict_data['map1_time'] * 100)
    print("Download time: " + f2(total_time) + "s")
    print("Download time percentage: " + perc + "%")
    total_time = float(dict_data['map1_upload_index']) + float(dict_data['map1_upload_map'])
    perc = f2(total_time / dict_data['map1_time'] * 100)
    print("Upload time: " + str(f2(total_time)) + "s")
    print("Upload time percentage: " + perc + "%")
    print("")
    
    print("*INDEX CORRECTION*")
    print("Average time: " + f2(dict_data['index_correction_time']) + "s")
    print("Download time: " + f2(dict_data['download_indexes']) + "s")
    perc = float(dict_data['download_indexes']) / float(dict_data['index_correction_time']) * 100
    print("Download time percentage: " + f2(perc) + "%")
    print("Upload time: " + f2(dict_data['upload_corrected_index']) + "s")
    perc = float(dict_data['upload_corrected_index']) / float(dict_data['index_correction_time']) * 100
    print("Upload time percentage: " + f2(perc) + "%")
    print("")
    
    print("*MAP STAGE 2*")
    print("Average time: " + f2(dict_data['map2_time']) + "s")
    total_time = float(dict_data['map2_download_fasta_chunk']) + float(dict_data['map2_download_map']) + float(dict_data['map2_download_corrected_index'])
    perc = total_time / dict_data['map2_time'] * 100
    print("Download time: " + f2(total_time) + "s")
    print("Download time percentage: " + f2(perc) + "%")
    print("Upload time: " + f2(dict_data['map2_upload_mpileup']))
    perc = float(dict_data['map2_upload_mpileup']) / float(dict_data['map2_time']) * 100
    print("Upload time percentage: " + f2(perc) + "%")
    print("")


def display_reduce(dict_data):
    print("**REDUCE PHASE**")
    print("*DISTRIBUTE INDEXES*")
    print("Average time: " + f2(dict_data['distribute_indexes_time']) + "s")
    print("S3 Select time: " + f2(dict_data['distribute_indexes_select']) + "s")
    perc = float(dict_data['distribute_indexes_select']) / dict_data['distribute_indexes_time'] * 100
    print("S3 Select time percentage: " + f2(perc) + "%")
    print("")
    
    print("*REDUCE*")
    print("Average time: " + f2(dict_data['reduce_function_time']) + "s")
    print("S3 SELECT time: " + f2(dict_data['reduce_func_select']) + "s")
    perc = float(dict_data['reduce_func_select']) / dict_data['reduce_function_time'] * 100
    print("S3 SELECT time percentage: " + f2(perc) + "%")
    print("Upload time: " + f2(dict_data['reduce_func_upload']) + "s")
    perc = float(dict_data['reduce_func_upload']) / dict_data['reduce_function_time'] * 100
    print("Upload time percentage: " + f2(perc) + "%")
    print("")
    
    print("*FINAL MERGE*")
    print("Average time: " + f2(dict_data['final_merge_total']) + "s")
    print("Download time: " + f2(dict_data['final_merge_download']) + "s")
    perc = float(dict_data['final_merge_download']) / dict_data['final_merge_total'] * 100
    print("Download time percentage: " + f2(perc) + "%")
    print("Upload time: " + f2(dict_data['final_merge_upload']) + "s")
    perc = float(dict_data['final_merge_upload']) / dict_data['final_merge_total'] * 100
    print("Upload time percentage: " + f2(perc) + "%")
    print("")
    
    
def display_all(dict_data):
    display_general(dict_data)
    display_map(dict_data)
    display_reduce(dict_data)
 

def plot():
    with open("stats/logs_stats.json") as json_read:
        data: dict = json.load(json_read)
    
    dict_data={}
    dict_data = fetch_general_data(data, dict_data)
    dict_data = fetch_map_data(data, dict_data)
    dict_data = fetch_reduce_data(data, dict_data)
    
    display_all(dict_data)


if __name__ == '__main__':
    plot()