import time
from lithops import Storage
import lithops
import subprocess as sp
from alignment_mapper import AlignmentMapper
from varcall_arguments import Arguments

def index_correction_map(setname: str, bucket: str, storage: Storage) -> int:
    """
    Corrects the index after the first map iteration. 
    All the set files must have the prefix "map_index_files/".
    Corrected indices will be stored with the prefix "correctedIndex/".
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
    storage.upload_file('/tmp/'+setname+'.txt', bucket, 'correctedIndex/'+setname+'.txt')
    return 0  

def delete_objects(storage: Storage, bucket: str, file_format: str):
    """
    Delete all objects with a given prefix inside a bucket
    """
    keys = storage.list_keys(bucket=bucket, prefix=file_format+"/")
    storage.delete_objects(bucket=bucket, key_list=keys)
        
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
        delete_objects(storage, args.bucket, args.file_format)

    print("Running Map Phase... " + str(len(iterdata)) + " functions")
    
    # Initizalize execution debug info
    start = time.time()
    
    if args.skip_map == "False":
        ###################################################################
        #### MAP: STAGE 1
        ###################################################################
        print("map phase - single iterdata set")
        fexec.map(map_func.map_alignment1, iterdata, timeout=int(args.func_timeout_map))
        first_map_results = fexec.get_result()
        
        ###################################################################
        #### MAP: GENERATE CORRECTED INDEXES
        ###################################################################
        # Generate the iterdata for index correction
        index_iterdata = []
        for i in range(num_chunks):
            index_iterdata.append({'setname': args.fq_seqname+'_fq'+str(i+1), 'bucket': str(args.bucket)})
        
        # Index correction
        fexec.map(index_correction_map, index_iterdata, timeout=int(args.func_timeout_map))
        corrections_results = fexec.get_result()
        
        ###################################################################
        #### MAP: STAGE 2
        ###################################################################
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
        fexec.map(map_func.map_alignment2, newiterdata, timeout=int(args.func_timeout_map))
        map_results = fexec.get_result()
        fexec.plot()
                            
    else:   # Skip map and get keys from previous run
        print("skipping map phase and retrieving existing keys")
        map_results = storage.list_keys(args.bucket, prefix="csv/")

    #End of map
    end = time.time()
    map_time = end - start

    #Delete intermediate files
    keys = storage.list_keys(args.bucket, "map_index_files/")
    for key in keys:
        storage.delete_object(args.bucket, key)
    keys = storage.list_keys(args.bucket, "correctedIndex/")
    for key in keys:
        storage.delete_object(args.bucket, key)
    keys = storage.list_keys(args.bucket, "filtered_map_files/")
    for key in keys:
        storage.delete_object(args.bucket, key)
    
    return map_time

    # TODO ADD REDUCE