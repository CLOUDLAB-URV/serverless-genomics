import time
from lithops import Storage
import lithops
import subprocess as sp
from alignment_mapper import AlignmentMapper
from varcall_arguments import Arguments
import os
import aux_functions as af

map1_cachefile = '/tmp/lithops_map1_checkpoint'
correction_cachefile = '/tmp/lithops_correction_checkpoint'
map2_cachefile = '/tmp/lithops_map2_checkpoint'

def index_correction_map(setname: str, bucket: str, storage: Storage) -> int:
    """
    Corrects the index after the first map iteration. 
    All the set files must have the prefix "map_index_files/".
    Corrected indices will be stored with the prefix "corrected_index/".
    """
    # Download all files related to this set
    filelist = storage.list_keys(bucket, "map_index_files/"+setname)
    for file in filelist:
        local_file = file.split("/")[-1]
        storage.download_file(bucket, file, '/tmp/'+local_file)
        
    # Execute correction scripts
    cmd = f'/function/bin/binary_reducer.sh /function/bin/merge_gem_alignment_metrics.sh 4 /tmp/{setname}* > /tmp/{setname}.intermediate.txt'
    result1 = sp.run(cmd, shell=True, check=True, universal_newlines=True)
    cmd2 = f'/function/bin/filter_merged_index.sh /tmp/{setname}.intermediate.txt /tmp/{setname}'
    result2 = sp.run(cmd2, shell=True, check=True, universal_newlines=True)
    
    # Upload corrected index to storage
    storage.upload_file('/tmp/'+setname+'.txt', bucket, 'corrected_index/'+setname+'.txt')
    return 0  
        
def map_reduce(args: Arguments, iterdata: list, map_func: AlignmentMapper, num_chunks: int) -> float:
    ###################################################################
    #### START OF MAP/REDUCE
    ###################################################################
    log_level = "DEBUG"
    
    # Initizalize storage and backend instances
    storage = Storage()
    s3 = storage.storage_handler.s3_client
    fexec = lithops.FunctionExecutor(log_level=log_level, runtime=args.runtime_id, runtime_memory=args.runtime_mem)

    if args.skip_map == "False":
        # Delete old files
        print("Deleting previous mapper outputs...")
        af.delete_objects(storage, args.bucket, args.file_format)

    print("Running Map Phase... " + str(len(iterdata)) + " functions")
    
    # Initizalize execution debug info
    start = time.time()
    
    if args.skip_map == "False":
        ###################################################################
        #### MAP: STAGE 1
        ###################################################################
        print("PROCESSING MAP: STAGE 1")
        
        #Load futures from previous execution
        map1_futures = af.load_cache(map1_cachefile)
        
        #Execute first map if futures were not found
        if(not map1_futures):
            map1_futures = fexec.map(map_func.map_alignment1, iterdata, timeout=int(args.func_timeout_map))
        
        #Get results either from the old futures or the new execution
        first_map_results = fexec.get_result(fs=map1_futures)
        
        #Dump futures into file
        af.dump_cache(map1_cachefile, map1_futures)
        
        
        ###################################################################
        #### MAP: GENERATE CORRECTED INDEXES
        ###################################################################
        print("PROCESSING INDEX CORRECTION")
        
        #Load futures from previous execution  
        correction_futures = af.load_cache(correction_cachefile)
        
        #Execute correction if futures were not found 
        if(not correction_futures):
            # Generate the iterdata for index correction
            index_iterdata = []
            for i in range(num_chunks):
                index_iterdata.append({'setname': args.fq_seqname+'_fq'+str(i+1), 'bucket': str(args.bucket)})
            
            # Execute index correction
            correction_futures = fexec.map(index_correction_map, index_iterdata, timeout=int(args.func_timeout_map))
        
        #Get results either from the old futures or the new execution    
        corrections_results = fexec.get_result(fs=correction_futures)
        
        #Dump futures into file
        af.dump_cache(correction_cachefile, correction_futures)
        
        
        ###################################################################
        #### MAP: STAGE 2
        ###################################################################
        print("PROCESSING MAP: STAGE 2")
        
        #Load futures from previous execution 
        map2_futures = af.load_cache(map2_cachefile)
        
        #Execute correction if futures were not found    
        if(not map2_futures):
            # Generate new iterdata
            newiterdata = []
            for worker in first_map_results:
                newiterdata.append({
                    'fasta_chunk': worker[0],
                    'fastq_chunk': worker[1],
                    'corrected_map_index_file': worker[2].split("-")[0]+".txt",
                    'filtered_map_file': worker[3],
                    'base_name': worker[4],
                    'old_id': worker[5]
                })
            
            # Execute second stage of map
            map2_futures = fexec.map(map_func.map_alignment2, newiterdata, timeout=int(args.func_timeout_map))
        
        #Get results either from the old futures or the new execution     
        map_results = fexec.get_result(fs=map2_futures)
        
        #Dump futures into file
        af.dump_cache(map2_cachefile, map2_futures)
                            
    else:   # Skip map and get keys from previous run
        print("skipping map phase and retrieving existing keys")
        map_results = storage.list_keys(args.bucket, prefix="csv/")

    #End of map
    end = time.time()
    map_time = end - start

    #Delete intermediate files
    af.delete_intermediate_files(storage, args, ['map_index_files/', 'corrected_index/', 'filtered_map_files/'], [map1_cachefile, map2_cachefile, correction_cachefile])
    
    return map_time

    # TODO ADD REDUCE