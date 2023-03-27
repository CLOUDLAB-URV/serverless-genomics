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


def plot_data_transfers(data, average=False):
    #ALIGN READS
    align_reads = data['pipeline']['alignReads_phase']['align_reads']['phases']['gem_indexer_mapper']['function_details']
    
    align_reads_download = 0
    align_reads_upload = 0
    for func in align_reads:
        k = list(func)[0]
        func = func[k]['data_sizes']
        keys = list(func)
        align_reads_download += func[keys[0]] + func[keys[1]]
        align_reads_upload += func[keys[2]] + func[keys[3]]
    if average:
        align_reads_download /= len(align_reads)
        align_reads_upload /= len(align_reads)
        
    
    #INDEX CORRECTION
    index_correction = data['pipeline']['alignReads_phase']['align_reads']['phases']['index_correction']['function_details']
    
    index_correction_downloads = 0
    index_correction_uploads = 0
    for func in index_correction:
        k = list(func)[0]
        func = func[k]['data_sizes']
        keys = list(func)
        last_k = keys.pop()
        for elem in keys:
            index_correction_downloads += func[elem]
        index_correction_uploads += func[last_k]
    if average:
        index_correction_downloads /= len(index_correction)
        index_correction_uploads /= len(index_correction)
    
    
    #MAP PHASE 2
    map_two = data['pipeline']['alignReads_phase']['align_reads']['phases']['filter_index_to_mpileup']['function_details']
    
    map_two_downloads = 0
    map_two_uploads = 0
    for func in map_two:
        k = list(func)[0]
        func = func[k]['data_sizes']
        keys = list(func)
        last_k = keys.pop()
        for elem in keys:
            map_two_downloads += func[elem]
        map_two_uploads += func[last_k]
    if average:
        map_two_downloads /= len(map_two)
        map_two_uploads /= len(map_two)
    
    #DISTRIBUTE INDEXES
    dist_indexes = data['pipeline']['reduce_phase']['reduce']['phases']['distribute_indexes']['function_details']
    
    dist_indexes_downloads = 0
    for func in dist_indexes:
        k = list(func)[0]
        func = func[k]['data_sizes']
        keys = list(func)
        first_k = keys.pop(0)
        dist_indexes_downloads += func[first_k]
    if average:
        dist_indexes_downloads /= len(dist_indexes)
    
    #REDUCE
    reduce = data['pipeline']['reduce_phase']['reduce']['phases']['reduce_function']['function_details']
    
    reduce_download = 0
    reduce_upload = 0
    for func in reduce:
        k = list(func)[0]
        func = func[k]['data_sizes']
        keys = list(func)
        reduce_download += func[keys[0]]
        reduce_upload += func[keys[1]]
    if average:
        reduce_download /= len(reduce)
        reduce_upload /= len(reduce)
    
    #PLOTS
    categories = ['Map One', 'Index Correction', 'Map Two', 'Dist. Indexes', 'Reduce']
    downloads = [align_reads_download, index_correction_downloads, map_two_downloads, dist_indexes_downloads, reduce_download]
    uploads = [align_reads_upload, index_correction_uploads, map_two_uploads, 0, reduce_upload]

    # Positions of bars on y-axis
    y_pos = np.arange(len(categories))

    # Width of bars
    bar_width = 0.1

    # Plotting the bars
    fig, ax = plt.subplots()
    ax.bar(y_pos - bar_width/2, downloads, width=bar_width, align='center', color='red', label='Downloads')
    ax.bar(y_pos + bar_width/2, uploads, width=bar_width, align='center', color='blue', label='Uploads')

    # Adding labels, title, and legend
    ax.set_xticks(y_pos)
    ax.set_xticklabels(categories)
    ax.set_ylabel('Data in MB')
    if average:
        ax.set_title('Average Data Transfers for one function')
    else:
        ax.set_title('Total Data Transfers')
    #ax.invert_xaxis()
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)

    fig = ax.get_figure()
    if average:
        fig.savefig("./stats/data_transfers_average.png", bbox_inches='tight', dpi=400)
    else:
        fig.savefig("./stats/data_transfers_total.png", bbox_inches='tight', dpi=400)

    return downloads, uploads


def display_data_transfer_values(downloads, uploads, average=False):
    if average:
        print("Average Data Transfers for one function")
    else:
        print("Total Data Transfers")
    
    print("Map Phase One: {:.2f}MB downloaded, {:.2f}MB uploaded. Generated data is {:.2f}% of the downloaded data.".format(downloads[0], uploads[0], (uploads[0]/downloads[0]*100)))
    print("Index Correction: {:.2f}MB downloaded, {:.2f}MB uploaded. Generated data is {:.2f}% of the downloaded data.".format(downloads[1], uploads[1], (uploads[1]/downloads[1]*100)))
    print("Map Phase Two: {:.2f}MB downloaded, {:.2f}MB uploaded. Generated data is {:.2f}% of the downloaded data.".format(downloads[2], uploads[2], (uploads[2]/downloads[2]*100)))
    print("Distribute Indexes: {:.2f}MB downloaded. No data is uploaded".format(downloads[3]))
    print("Reduce Function: {:.2f}MB downloaded, {:.2f}MB uploaded. Generated data is {:.2f}% of the downloaded data.".format(downloads[4], uploads[4], (uploads[4]/downloads[4]*100)))


def plot_stages(data):
    #Fetch times
    fastq_preprocessing = data['pipeline']['preprocess_phase']['preprocess']['subprocesses_fastq']['get_data_frame_parquet']['execution_time'] \
        + data['pipeline']['preprocess_phase']['preprocess']['subprocesses_fastq']['prepare_fastq_chunks']['execution_time']
    fasta_preprocessing = data['pipeline']['preprocess_phase']['preprocess']['subprocesses_fasta']['prepare_fasta_chunks']['execution_time']
    
    map_one = data['pipeline']['alignReads_phase']['align_reads']['phases']['gem_indexer_mapper']['execution_time']
    index_correction = data['pipeline']['alignReads_phase']['align_reads']['phases']['index_correction']['execution_time']
    map_two = data['pipeline']['alignReads_phase']['align_reads']['phases']['filter_index_to_mpileup']['execution_time']
    
    dist_indexes = data['pipeline']['reduce_phase']['reduce']['phases']['distribute_indexes']['execution_time']
    reduce = data['pipeline']['reduce_phase']['reduce']['phases']['reduce_function']['execution_time']
    
    merge = data['pipeline']['reduce_phase']['reduce']['phases']['final_merge']['execution_time']
    
    #Merge into bars
    preprocessing_bar = [fastq_preprocessing, fasta_preprocessing, 0, 0, 0, 0, 0, 0]
    map_bar = [0, 0, map_one, index_correction, map_two, 0, 0, 0]
    reduce_bar = [0, 0, 0, 0, 0, dist_indexes, reduce, 0]
    merge_bar = [0, 0, 0, 0, 0, 0, 0, merge]
    
    x = ['Preprocessing', 'Map', 'Reduce', 'Merge']
    bars = [[0 for j in range(4)] for i in range(8)]
    i = 0
    for x1, x2, x3, x4 in zip(preprocessing_bar, map_bar, reduce_bar, merge_bar):
        bars[i][0] = x1
        bars[i][1] = x2
        bars[i][2] = x3
        bars[i][3] = x4
        i += 1
    
    values = np.array(bars)
    fig, ax = plt.subplots()
    
    colors = ['red', 'orange', 'forestgreen', 'limegreen', 'lime', 'navy', 'blue', 'magenta']
    
    # Stacked bar chart with loop
    for i in range(values.shape[0]):
        ax.bar(x, values[i], bottom = np.sum(values[:i], axis = 0), color=colors[i])
        
    ax.legend(["fastq_preprocessing","fasta_preprocessing","map_stage_one","index_correction",\
        "map_stage_two","distribute_indexes","reduce_function","merge"], \
        bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
    ax.set_title(f"{data['fasta_path'].split('/')[-1]} - 25.7MB - {data['fasta_chunks']} chunks")
    ax.set_xlabel("Runtime Mem: 2GB - Reducer Runtime Mem: 8GB")
    ax.set_ylabel("Time in seconds")
    
    fig = ax.get_figure()
    fig.suptitle(f"{data['fastq_path'].split('/')[-1]} - 685.3MB - {data['fastq_chunks']} chunks")
    fig.savefig("./stats/total_times_by_stage.png", bbox_inches='tight', dpi=400)
    
    return data['pipeline']['execution_time']


def plot_stages_num(data):
    #Fetch functions
    map_one = data['pipeline']['alignReads_phase']['align_reads']['phases']['gem_indexer_mapper']['function_details']
    index_correction = data['pipeline']['alignReads_phase']['align_reads']['phases']['index_correction']['function_details']
    map_two = data['pipeline']['alignReads_phase']['align_reads']['phases']['filter_index_to_mpileup']['function_details']
    
    dist_indexes = data['pipeline']['reduce_phase']['reduce']['phases']['distribute_indexes']['function_details']
    reduce = data['pipeline']['reduce_phase']['reduce']['phases']['reduce_function']['function_details']
    
    merge = data['pipeline']['reduce_phase']['reduce']['phases']['final_merge']['function_details']
    
    #Merge into bars
    map_bar = [len(map_one), len(index_correction), len(map_two), 0, 0, 0]
    reduce_bar = [0, 0, 0, len(dist_indexes), len(reduce), 0]
    merge_bar = [0, 0, 0, 0, 0, len(merge)]
    
    x = ['Map', 'Reduce', 'Merge']
    bars = [[0 for j in range(3)] for i in range(6)]
    i = 0
    for x1, x2, x3 in zip(map_bar, reduce_bar, merge_bar):
        bars[i][0] = x1
        bars[i][1] = x2
        bars[i][2] = x3
        i += 1
    
    values = np.array(bars)
    fig, ax = plt.subplots()
    
    colors = ['forestgreen', 'limegreen', 'lime', 'navy', 'blue', 'magenta']
    
    # Stacked bar chart with loop
    for i in range(values.shape[0]):
        ax.bar(x, values[i], bottom = np.sum(values[:i], axis = 0), color=colors[i])
        
    ax.legend(["map_stage_one","index_correction",\
        "map_stage_two","distribute_indexes","reduce_function","merge"], \
        bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
    ax.set_title(f"{data['fasta_path'].split('/')[-1]} - 25.7MB - {data['fasta_chunks']} chunks")
    ax.set_xlabel("Runtime Mem: 2GB - Reducer Runtime Mem: 8GB")
    ax.set_ylabel("Number of functions")
    
    fig = ax.get_figure()
    fig.suptitle(f"{data['fastq_path'].split('/')[-1]} - 685.3MB - {data['fastq_chunks']} chunks")
    fig.savefig("./stats/total_functions_by_stage.png", bbox_inches='tight', dpi=400)
    
    total = {
        'map_one': len(map_one),
        'index_correction': len(index_correction),
        'map_two': len(map_two),
        'dist_indexes': len(dist_indexes),
        'reduce': len(reduce),
        'merge': len(merge)
    }
    
    return total

def display_num_func(total):
    print(f'Map One: {total["map_one"]} functions launched.')
    print(f'Index Correction: {total["index_correction"]} functions launched.')
    print(f'Map Two: {total["map_two"]} functions launched.')
    print(f'Distribute Indexes: {total["dist_indexes"]} functions launched.')
    print(f'Reduce: {total["reduce"]} functions launched.')
    print(f'Merge: {total["merge"]} functions launched.')


if __name__ == '__main__':
    with open("/home/agabriel/Downloads/logs_stats.json") as json_read:
        data: dict = json.load(json_read)
    
    if not os.path.exists('stats'):
        os.makedirs('stats')
    
    t = plot_stages_num(data)
    display_num_func(t)