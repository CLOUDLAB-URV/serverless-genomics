#from codecs import iterdecode
from multiprocessing.sharedctypes import Value
from operator import methodcaller
from sys import prefix
import time
from lithops import Storage
import boto3
import lithops
import lithopsgenetics.auxiliaryfunctions as af
from PriceEstimator import PriceEstimator
from botocore.client import Config
import subprocess as sp


class MapReduce:

    def create_intermediate_keys(self, map_results):

        intermediate_keys = []
        aux = []

        # Create a list with all the intermediate keys -> [ [fasta_split1_keys], [fasta_split2_keys], ...]
        # Group keys generated in map phase into fasta_split sets
        # Start with first split number
        split_number = map_results[0].split("/")
        split_number = split_number[1].split("_")
        split_number = split_number[1].split(".")
        split_number = split_number[0]
        #af.printl("split number for intermediate keys: "+str(split_number), self.stage, self.id, self.debug)


        for key in map_results:
            #af.printl("key from map function : "+key, self.stage, self.id, self.debug)
            #af.printl("split number in for loop: "+str(split_number), self.stage, self.id, self.debug)
            split = key.split("/")
            split = split[1].split("_")
            split = split[1].split(".")
            # if fasta_split in map key corresponds to current fasta_split no., add to set
            if(split[0] == split_number):
                aux.append(key)
            # if fasta split in map key does not correspond, save existing set to intermediate keys, and start a new set with the new split number from the map key.
            else:
                intermediate_keys.append(aux)
                aux = []
                aux.append(key)
                split_number = split[0]
                

        intermediate_keys.append(aux)
        return intermediate_keys
        

    #Create index (mpileup position column) dictionary -> {index : #rows}
    def get_s3_select_indexes(self, bucket, keys, file_format):
        s3 = boto3.client('s3')

        count_indexes = {}

        if file_format == "csv":
            expression = "SELECT cast(s._2 as int) FROM s3object s"
            input_serialization = {'CSV': {}, 'CompressionType': 'NONE'}

        elif file_format == "parquet":
            expression = "SELECT s.\"1\" FROM s3object s"
            input_serialization = {'Parquet': {}, 'CompressionType': 'NONE'}
        else:
            return "ERROR: Invalid format"
        
        # 1. RETRIEVE MPILEUP POSITIONS WITH S3 SELECT
        start=time.time()
        if isinstance(keys, str):
            keys = [keys]
        for key in keys:
            resp = s3.select_object_content(
                Bucket=bucket,
                Key=key,
                ExpressionType='SQL',
                Expression=expression,
                InputSerialization = input_serialization,
                OutputSerialization = {'CSV': {}}
            )

            data = ""
            for event in resp['Payload']:
                if 'Records' in event:
                    records = event['Records']['Payload'].decode("UTF-8")
                    data = data + records

            str_indexes = data.split("\n")
            str_indexes.pop()
            int_indexes = list(map(int, str_indexes))

            for index in int_indexes:
                if index in count_indexes:
                    count_indexes[index] += 1
                else:
                    count_indexes[index] = 1
        af.printl("s3 select: execution_time: "+str(time.time()-start)+": s",self.stage, self.id, self.debug)
        print("s3 select2: execution_time: "+str(time.time()-start)+": s")
        # 2. COUNT INDEXES
        start=time.time()
        count_indexes = dict(sorted(count_indexes.items()))
        af.printl("count indexes: execution_time: "+str(time.time()-start)+": ns",self.stage, self.id, self.debug)
        
        indexes = 0
        #MAX rows a reducer can process
        max_index = 20000000
        index_list = []

        # 3. GENERATE RANGES
        start=time.time()
        # summing all index/position counts while the count is lower than the maximum index
        for key in count_indexes:
            if indexes + count_indexes[key] < max_index:
                indexes = indexes + count_indexes[key]
                index = key
            else: # append the last index below max_index as end value in range, and start a new range.
                indexes = 0
                index_list.append(index)

        last = list(count_indexes)[-1]
        index_list.append(last)
        af.printl("generating ranges: execution_time: "+str(time.time()-start)+": ns",self.stage, self.id, self.debug)
        
        return index_list 

    def calculate_ranges(self, bucket, rows_keys):

        storage = Storage()

        index_ranges = {}
        row_key_count=0
        for key in rows_keys:
            row_key_count+=1
            print("row key no. " +str(row_key_count) +": " + str(key))
            rows = storage.get_object(bucket=bucket, key=key).decode("UTF-8")
            rows = rows.split("\n")
            rows.pop(-1)
            print("number of rows: "+ str(len(rows)))
            print("row 1: " + str(rows[1]))
            print("row 1 index : " + str(rows[1].split(":")[0]) )

            for row in rows:
                index = row.split(":")
                if index[0] in index_ranges:
                    index_ranges[int(index[0])] += int(index[1])
                else:
                    index_ranges[int(index[0])] = int(index[1])
        
        index_ranges = dict(sorted(index_ranges.items()))

        range_index = []
        cur_rows = 0
        aux_range = 0
        for index, rows in index_ranges.items():
            cur_rows = cur_rows + int(rows)
            if (cur_rows > 20000000):
                range_index.append(aux_range)
                cur_rows = 0
            
            aux_range = index

        range_index.append(aux_range)

        return range_index

    def create_multipart(self, s3, bucket, key):
        mpu = s3.create_multipart_upload(
            Bucket=bucket,
            Key=key
        )
        return mpu['UploadId']


    def complete_multipart(self, s3, bucket, keys, mpu_ids, parts):
        results = []

        for key, mpu_id in zip(keys, mpu_ids):
            mpu_part = []
            remove = 0

            for part in parts:
                if mpu_id == part['mpu_id']:
                    mpu_part.append({"PartNumber" : part["PartNumber"], "ETag" : part["ETag"]})
                    remove = remove + 1
                else:
                    break

            result = s3.complete_multipart_upload(
                Bucket = bucket,
                Key = key,
                UploadId = mpu_id,
                MultipartUpload = {"Parts": mpu_part}
                )
            results.append(result)
    
            for _ in range(remove):
                parts.pop(0)
                
        return results

    #Create the iterdata for the reduce function
    def create_iterdata_reducer(self, mpileup_keys, indexes_ranges, mpu_ids, mpu_keys, file_format, buffer_size):
        iterdata = []

        for keys, indexes, mpu_id, mpu_key in zip(mpileup_keys, indexes_ranges, mpu_ids, mpu_keys):
            start = 1
            n_part = 1

            for index in indexes:
                data = {
                    "key" : keys, 
                    "range" : { "start" : start,
                                "end" : int(index)
                            },
                    "mpu_id" : mpu_id,
                    "n_part" : n_part,
                    "mpu_key" : mpu_key,
                    "file_format" : file_format,
                    "buffer_size" : buffer_size
                }
                iterdata.append(data)
                n_part += 1
                start = int(index) + 1

        return iterdata

    def final_merge(self, storage, bucket, mpu_id, mpu_key, key, n_part):

        sinple_out = storage.get_object(bucket=bucket, key=key)

        #Upload part
        s3 = boto3.client('s3')
        part = s3.upload_part(
            Body = sinple_out,
            Bucket = bucket,
            Key = mpu_key,
            UploadId = mpu_id,
            PartNumber = n_part
        )
        
        return {"PartNumber" : n_part, "ETag" : part["ETag"], "mpu_id": mpu_id}


    def delete_objects(self, storage, bucket, file_format):
        keys = storage.list_keys(bucket=bucket, prefix=file_format+"/")
        storage.delete_objects(bucket=bucket, key_list=keys)


    def index_correction_map(id, setname, bucket, storage):  
        filelist = storage.list_keys(bucket, "map_index_files/"+setname)
        for file in filelist:
            local_file = file.split("/")[-1]
            storage.download_file(bucket, file, '/tmp/'+local_file)
        
        cmd = f'/function/bin/binary_reducer.sh /function/bin/merge_gem_alignment_metrics.sh 4 /tmp/{setname}* > /tmp/{setname}.intermediate.txt'
        result1 = sp.run(cmd, shell=True, check=True, universal_newlines=True)
        cmd2 = f'/function/bin/filter_merged_index.sh /tmp/{setname}.intermediate.txt /tmp/{setname}'
        result2 = sp.run(cmd2, shell=True, check=True, universal_newlines=True)
        
        storage.upload_file('/tmp/'+setname+'.txt', bucket, 'correctedIndex/'+setname+'.txt')
        return 0
    
    
    def __init__(self, map_func, map_func2, reduce_func, create_sinple_key, runtime, runtime_memory, runtime_memory_r, buffer_size, file_format, func_timeout_map, func_timeout_reduce,log_level, bucket, stage, id, debug, skip_map, workers, method, fastq_set_n, run_id, iterdata_n, num_chunks, fq_seqname):
        self.map_func = map_func
        self.map_func2 = map_func2
        self.reduce_func = reduce_func
        self.create_sinple_key = create_sinple_key
        self.runtime = runtime
        self.runtime_memory = runtime_memory
        self.runtime_memory_r = runtime_memory_r
        self.buffer_size = buffer_size
        self.file_format = file_format
        self.func_timeout_map = func_timeout_map
        self.func_timeout_reduce = func_timeout_reduce
        self.log_level = log_level
        self.bucket = bucket
        self.stage = stage
        self.id = id
        self.debug = debug
        self.skip_map = skip_map
        self.method = method
        self.workers = workers
        self.fastq_set_n = fastq_set_n
        self.run_id = run_id
        self.iterdata_n = iterdata_n
        self.num_chunks = num_chunks
        self.fq_seqname = fq_seqname
    
    
    

    def __call__(self, iterdata, iterdata_sets):

        storage = Storage()
        s3 = storage.storage_handler.s3_client

        # log_level=self.log_level
        fexec = lithops.FunctionExecutor(log_level=self.log_level, runtime=self.runtime, runtime_memory=self.runtime_memory)

        if self.skip_map == "False":
            print("Deleting previous mapper outputs...")
            self.delete_objects(storage, self.bucket, self.file_format)

        print("Running Map Phase... " + str(len(iterdata)) + " functions")
        
        stage="A"
        id="X"
        start = time.time()
        map_results=[]
        if self.skip_map == "False":
            if not iterdata_sets:
                print("map phase - single iterdata set")
                
                # Start first part of map
                fexec.map(self.map_func, iterdata, timeout=2400)
                first_map_results = fexec.get_result()
                
                # Generate index correction iterdata
                index_iterdata = []
                for i in range(self.num_chunks):
                    index_iterdata.append({'setname': self.fq_seqname+'_fq'+str(i+1), 'bucket': str(self.bucket)})
                
                # Index correction
                fexec.map(self.index_correction_map, index_iterdata)
                corrections_results = fexec.get_result()
                
                # Generate new iterdata
                newiterdata = []
                for worker in first_map_results:
                    newiterdata.append({
                        'fasta_chunk': worker[0],
                        'fastq_chunk': worker[1],
                        'corrected_map_index_file': worker[2].split("-")[0]+".txt",
                        'filtered_map_file': worker[3],
                        'base_name': worker[4],
                        'gem_index_present': worker[5],
                        'old_id': worker[6]
                    })
                
                # Execute second part of map
                fexec.map(self.map_func2, newiterdata, timeout=2400)
                
                # Results
                map_results = fexec.get_result()
                fexec.plot()
            else:
                print("map phase - iterdata with "+str(len(iterdata_sets))+" sets")
                count=0
                for iterdata in iterdata_sets:
                    count+=1
                    #print("map phase - iterdata set no. "+str(count)+" with "+str(len(iterdata))+" elements")
                    #print("\nprinting iterdata")
                    #print(str(iterdata))
                    #print("\nprinting iterdata elements")
                    #for el in iterdata:
                        #print(el)
                    fexec.map(self.map_func, iterdata, timeout=2400)
                    print("map phase - finished iterdata set no. "+str(count)+" with "+str(len(iterdata))+" elements")
                    print("getting results")
                    map_results_part = fexec.get_result()
                    fexec.plot()
                    print("appending partial results to map results")
                    for el in map_results_part:
                        map_results.append(el)
                    
        else:
            print("skipping map phase and retrieving existing keys")
            map_results = storage.list_keys(self.bucket, prefix="csv/")
            
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

        end = time.time()
        map_time = end - start
        print(f'{stage}:{id}: map phase: execution_time: {map_time}: s')
        pEst = PriceEstimator(fexec) 
        runtime_mem_gb= int(self.runtime_memory)/1024
        print('map: lambda: cost:'+ str(pEst.lambda_calc(runtime_mem_gb)))
        print('map: s3: cost:'+ str(pEst.s3_calc_multiple_fexec(20)))
        print("Map Phase Finished...")

        return map_time, 0, 0, 0 # Skip reduce for testing purposes

        #print("\nmap results\n"+str(map_results))

        #--------------------------------------------------------------------------
        # REDUCE STAGE 1
        print("Reduce stage 1: Creating Intermediate Keys and mpileup index ranges...")
        start = time.time()
        stage = "B"
        id = "X"

        object_keys = []
        rows_keys = []
        
        for results in map_results:
            object_keys.append(results[0])
            rows_keys.append(results[1])
        del map_results

        intermediate_keys = self.create_intermediate_keys(object_keys)
        rows_keys = self.create_intermediate_keys(rows_keys)
        
        af.printl("Expected reduce functions - intermediate keys (STAGE 1): " + str(len(intermediate_keys)), stage, id, self.debug)

        sinple_key = self.create_sinple_key()
        final_sinple_key = []
        final_sinple_key.append(sinple_key.pop(-1))
        
        sinple_key_count= 0
        for x in sinple_key:
            sinple_key_count+=1
            #print("sinple key no. "+str(sinple_key_count) + ": " + str(x))
        
        #MultipartUpload Ids -> [id, id, ...]
        mpu_ids = []
        for key in sinple_key:
            mpu_ids.append(self.create_multipart(s3, self.bucket, key))

        ############
        if(self.method == "select"):
            indexes = []
            for keys in object_keys:
                data = {
                    "bucket": self.bucket ,
                    "keys": keys,
                    "file_format" : self.file_format
                }
                indexes.append(data)

            fexec = lithops.FunctionExecutor(log_level=self.log_level, runtime=self.runtime)
            fexec.map(self.get_s3_select_indexes, indexes, timeout=int(self.func_timeout_reduce))
            results = fexec.get_result()
            pEst = PriceEstimator(fexec) 
            runtime_mem_gb= 8192/1024
            print('reduce 1 (select): lambda: cost:'+ str(pEst.lambda_calc(runtime_mem_gb)))
            print('reduce 1 (select): s3: cost:' +str(pEst.s3_calc_multiple_fexec(2)))
            fexec.plot()

        ############
        elif(self.method == "manual"):

            indexes = []
            for keys in rows_keys:
                data = {
                    "bucket": self.bucket,
                    "rows_keys": keys,
                }
                indexes.append(data)
            
            fexec = lithops.FunctionExecutor(log_level=self.log_level, runtime=self.runtime)
            fexec.map(self.calculate_ranges, indexes, timeout=int(self.func_timeout_reduce))
            fexec.wait()
            results = fexec.get_result()

            pEst = PriceEstimator(fexec)
            runtime_mem_gb= 4096/1024
            print('reduce 1 (manual): lambda: cost:'+ str(pEst.lambda_calc(runtime_mem_gb)))
            print('reduce 1 (manual): s3: cost:' + str(pEst.s3_calc_multiple_fexec(0)))

            fexec.plot()
        #-----------        #Iterdata for the reducers
        iterdata = self.create_iterdata_reducer(intermediate_keys, results, mpu_ids, sinple_key, self.file_format, self.buffer_size)
        
        end = time.time()
        creating_keys_time = end - start
        print(f'{stage}:{id}: reduce 1 phase: execution_time: {creating_keys_time}: s')
        #for i in iterdata:
        #    for x in i.items():
                #print(x)

        
        
        #--------------------------------------------------------------------------

        
        print("Reduce Phase 2: generating sinple output from mpileup")
        start = time.time()
        stage = "C"
        af.printl("Expected reduce functions (STAGE 2): " + str(len(iterdata)), stage, id, self.debug)

        fexec = lithops.FunctionExecutor(log_level=self.log_level, runtime=self.runtime)
        fexec.map(self.reduce_func, iterdata, timeout=int(self.func_timeout_reduce))
        results = fexec.get_result()
        fexec.plot()
        pEst = PriceEstimator(fexec) 
    
        runtime_mem_gb= int(self.runtime_memory_r)/1024
        print('reduce 2: lambda: cost:'+ str(pEst.lambda_calc(runtime_mem_gb)))
        print('reduce 2: s3: cost:' +  str(pEst.s3_calc_multiple_fexec(4)))
        end = time.time()
        reduce_time = end - start
        print(f'{stage}:{id}: reduce 2 phase: execution_time: {reduce_time}: s')

        #--------------------------------------------------------------------------
        print("reduce 3: multipart upload")
        stage = "D"
        #results = self.complete_multipart(s3, self.bucket, sinple_key, mpu_ids, results)
        start = time.time()

        n_parts = len(sinple_key)
        mpu_id = self.create_multipart(s3, self.bucket, final_sinple_key[0])

        print("creating multipart iterdata")
        part = 1
        iterdata = []
        while part <= n_parts:
            # print("bucket: " + str(self.bucket))
            # print("mpu_id: " + str(mpu_id))
            # print("mpu_key: " + final_sinple_key[0])
            # print("key: " + str(sinple_key[part-1]))
            # print("n_part: " + str(part))
            data = {
                "bucket": self.bucket,
                "mpu_id": mpu_id,
                "mpu_key": final_sinple_key[0],
                "key": sinple_key[part-1],
                "n_part": part
            }
            iterdata.append(data)
            part += 1


        fexec = lithops.FunctionExecutor(log_level=self.log_level, runtime=self.runtime)
        fexec.map(self.final_merge, iterdata, timeout=int(self.func_timeout_reduce))
        results = fexec.get_result()
        fexec.plot()
        pEst = PriceEstimator(fexec) 

        mpu_ids = []
        mpu_ids.append(mpu_id)
        results = self.complete_multipart(s3, self.bucket, final_sinple_key, mpu_ids, results)



        for x in results:
            print(x)

        fexec.plot()
        runtime_mem_gb= int(self.runtime_memory_r)/1024
        print('output merge: lambda: cost:'+ str(pEst.lambda_calc(runtime_mem_gb)))
        print('output merge: s3: cost:' +  str(pEst.s3_calc_multiple_fexec(1)))

        end = time.time()
        multipart_upload_time = end - start
        print(f'{stage}:{id}: reduce 3 phase: execution_time: {reduce_time}: s')

        print("Terminating MapReducer...\n")
        return map_time, creating_keys_time, reduce_time, multipart_upload_time
        
    