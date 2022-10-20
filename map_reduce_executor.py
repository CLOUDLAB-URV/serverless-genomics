import time
from lithops import Storage
import lithops
import subprocess as sp
from map_functions import MapFunctions

class MapReduce:
    """
    The objective of this class is to coordinate all the steps involving the Processing stage of the pipeline. This
    includes executing the two maps, the reduce and the index correction functions.
    """
    def __init__(self, map_func, runtime, runtime_memory, runtime_memory_r, buffer_size, file_format, func_timeout_map, 
                 func_timeout_reduce,log_level, bucket, stage, id, debug, skip_map, method, iterdata_n, num_chunks, fq_seqname):
        self.map_func: MapFunctions = map_func
        self.runtime = runtime
        self.runtime_memory = runtime_memory
        self.runtime_memory_r = runtime_memory_r
        self.buffer_size = buffer_size
        self.file_format = file_format
        self.func_timeout_map = int(func_timeout_map)
        self.func_timeout_reduce = int(func_timeout_reduce)
        self.log_level = log_level
        self.bucket = bucket
        self.stage = stage
        self.id = id
        self.debug = debug
        self.skip_map = skip_map
        self.method = method
        self.iterdata_n = iterdata_n
        self.num_chunks = num_chunks
        self.fq_seqname = fq_seqname


    def delete_objects(self, storage, bucket, file_format):
        """
        Delete all objects with a given prefix inside a bucket
        """
        keys = storage.list_keys(bucket=bucket, prefix=file_format+"/")
        storage.delete_objects(bucket=bucket, key_list=keys)


    def index_correction_map(id, setname, bucket, storage):
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


    def __call__(self, iterdata):
        ###################################################################
        #### START OF MAP/REDUCE
        ###################################################################
        
        # Initizalize storage and backend instances
        storage = Storage()
        s3 = storage.storage_handler.s3_client
        fexec = lithops.FunctionExecutor(log_level=self.log_level, runtime=self.runtime, runtime_memory=self.runtime_memory)

        if self.skip_map == "False":
            # Delete old files
            print("Deleting previous mapper outputs...")
            self.delete_objects(storage, self.bucket, self.file_format)

        print("Running Map Phase... " + str(len(iterdata)) + " functions")
        
        # Initizalize execution debug info
        stage="A"
        id="X"
        start = time.time()
        map_results=[]
        
        if self.skip_map == "False":
            ###################################################################
            #### MAP: STAGE 1
            ###################################################################
            print("map phase - single iterdata set")
            fexec.map(self.map_func.map_alignment1, iterdata, timeout=self.func_timeout_map)
            first_map_results = fexec.get_result()
            
            ###################################################################
            #### MAP: GENERATE CORRECTED INDEXES
            ###################################################################
            # Generate the iterdata for index correction
            index_iterdata = []
            for i in range(self.num_chunks):
                index_iterdata.append({'setname': self.fq_seqname+'_fq'+str(i+1), 'bucket': str(self.bucket)})
            
            # Index correction
            fexec.map(self.index_correction_map, index_iterdata, timeout=self.func_timeout_map)
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
            fexec.map(self.map_func.map_alignment2, newiterdata, timeout=self.func_timeout_map)
            map_results = fexec.get_result()
            fexec.plot()
                             
        else:   # Skip map and get keys from previous run
            print("skipping map phase and retrieving existing keys")
            map_results = storage.list_keys(self.bucket, prefix="csv/")

        #End of map
        end = time.time()
        map_time = end - start

        #Delete intermediate files
        keys = storage.list_keys(self.bucket, "map_index_files/")
        for key in keys:
            storage.delete_object(self.bucket, key)
        keys = storage.list_keys(self.bucket, "correctedIndex/")
        for key in keys:
            storage.delete_object(self.bucket, key)
        keys = storage.list_keys(self.bucket, "filtered_map_files/")
        for key in keys:
            storage.delete_object(self.bucket, key)
        
        return map_time 
    
        # TODO ADD REDUCE