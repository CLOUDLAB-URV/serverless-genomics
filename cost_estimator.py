import json
import os
import csv
import lithops

def cost_estimation(data, cost_mem, cost_reduce_mem, select_scan_cost, select_return_cost, bucket, storage):
    #Map One
    map_one_sum = 0
    
    function_details = data['pipeline']['alignReads_phase']['align_reads']['phases']['gem_indexer_mapper']['function_details']
    for elem in function_details:
        k = list(elem.keys())[0]
        timestamps = elem[k]['timestamps']
        map_one_sum += (timestamps['end'] - timestamps['start'])
        
    #Index Correction
    index_sum = 0
    
    function_details = data['pipeline']['alignReads_phase']['align_reads']['phases']['index_correction']['function_details']
    for elem in function_details:
        k = list(elem.keys())[0]
        timestamps = elem[k]['timestamps']
        index_sum += (timestamps['end'] - timestamps['start'])
    
    #Map Two
    map_two_sum = 0
    
    function_details = data['pipeline']['alignReads_phase']['align_reads']['phases']['filter_index_to_mpileup']['function_details']
    for elem in function_details:
        k = list(elem.keys())[0]
        timestamps = elem[k]['timestamps']
        map_two_sum += (timestamps['end'] - timestamps['start'])
        
    #Distribute Indexes
    dist_indexes_sum = 0
    
    function_details = data['pipeline']['reduce_phase']['reduce']['phases']['distribute_indexes']['function_details']
    for elem in function_details:
        k = list(elem.keys())[0]
        timestamps = elem[k]['timestamps']
        dist_indexes_sum += (timestamps['end'] - timestamps['start'])
        
    #Reduce
    reduce_sum = 0
    
    function_details = data['pipeline']['reduce_phase']['reduce']['phases']['reduce_function']['function_details']
    for elem in function_details:
        k = list(elem.keys())[0]
        timestamps = elem[k]['timestamps']
        reduce_sum += (timestamps['end'] - timestamps['start'])
        
    #Merge
    merge_sum = 0
    
    function_details = data['pipeline']['reduce_phase']['reduce']['phases']['final_merge']['function_details']
    for elem in function_details:
        for k in elem:
            merge_sum += elem[k]['execution_time']
        
    #S3 Select Operations
    dist_select_scan = 0
    dist_select_returned = 0
    
    function_details = data['pipeline']['reduce_phase']['reduce']['phases']['distribute_indexes']['function_details']
    for elem in function_details:
        k = list(elem.keys())[0]
        dist_select_returned += int(elem[k]['data_sizes']['total_data_from_select']) / 1000
        for file_k in elem[k]['data_sizes']['keys']:
            resp = storage.head_object(bucket=bucket, key=file_k)
            dist_select_scan += (int(resp['content-length']) / (1000*1000*1000))
    
    reduce_select_scan = 0
    reduce_select_returned = 0
    
    function_details = data['pipeline']['reduce_phase']['reduce']['phases']['reduce_function']['function_details']
    for elem in function_details:
        k = list(elem.keys())[0]
        k2 = list(elem[k]['data_sizes'].keys())[0]
        reduce_select_returned += int(elem[k]['data_sizes'][k2]) / 1000
        for file_k in elem[k]['data_sizes']['keys']:
            resp = storage.head_object(bucket=bucket, key=file_k)
            reduce_select_scan += (int(resp['content-length']) / (1000*1000*1000))

    
    # Open a new CSV file in write mode
    with open('stats/costs.csv', 'w', newline='') as file:
        writer = csv.writer(file)

        # Write the header row
        writer.writerow(['Stage', 'Cost'])

        # Write each stage and its corresponding cost
        writer.writerow(['map_one', map_one_sum*cost_mem])
        writer.writerow(['index_correction', index_sum*cost_mem])
        writer.writerow(['map_two', map_two_sum*cost_mem])
        
        writer.writerow(['dist_indexes', dist_indexes_sum*cost_reduce_mem])
        writer.writerow(['reduce', reduce_sum*cost_reduce_mem])
        writer.writerow(['merge', merge_sum*cost_reduce_mem])
        
        writer.writerow(['dist_indexes_select', dist_select_scan*select_scan_cost+dist_select_returned*select_return_cost])
        writer.writerow(['reduce_indexes_select', reduce_select_scan*select_scan_cost+reduce_select_returned*select_return_cost])
        
        total = map_one_sum*cost_mem + index_sum*cost_mem + map_two_sum*cost_mem \
            + dist_indexes_sum*cost_reduce_mem + reduce_sum*cost_reduce_mem + merge_sum*cost_reduce_mem \
            + dist_select_scan*select_scan_cost+dist_select_returned*select_return_cost + reduce_select_scan*select_scan_cost+reduce_select_returned*select_return_cost
        writer.writerow(['total', total])
    

if __name__ == '__main__':
    with open("/home/agabriel/Downloads/logs_stats.json") as json_read:
        data: dict = json.load(json_read)
    
    if not os.path.exists('stats'):
        os.makedirs('stats')
        
    storage = lithops.Storage()

    cost_estimation(data, 0.0000333, 0.0001333, 0.002, 0.0007, 'agabriel-data', storage)