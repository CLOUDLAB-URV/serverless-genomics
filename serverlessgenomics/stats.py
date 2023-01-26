from lithops import Storage
import time
import json
from copy import deepcopy

class Stats:
    def __init__(self):
        self.__tmp_registrer={}
        self.__stats={}

    def timer_start(self, script, extra_time=None):
        if script in self.__tmp_registrer:
            print(f'WARNING: the counter of the timer \"{script}\" was already running, it will restart.')
        elif script in self.__stats and "execution_time" in self.__stats[script]:
            raise Exception(f'The timer of \"{script}\" already existed, choose another name or remove the existing timer first.')
 
        self.__tmp_registrer[script] = time.perf_counter() if extra_time is None else extra_time + time.perf_counter()
    
    def timer_stop(self, script):
        end_time = time.perf_counter()
        
        if script in self.__tmp_registrer:
            if script not in self.__stats:
                self.__stats[script] = {} 
            self.__stats[script]["execution_time"] = end_time - self.__tmp_registrer[script]
            del self.__tmp_registrer[script]
        else:
            raise Exception(f'You can not execute this function before "timer_start".')
    
    def store_size_data(self, name_data, size, script=None):
        if script is None: 
                dictionary = self.__stats
        else:                
            if script not in self.__stats:
                self.__stats[script] = {} 
            dictionary = self.__stats[script]

        if name_data in dictionary:
            raise Exception(f'The key \"{name_data}\" contains data, choose another name or remove the key first.')
        else:
            dictionary[name_data] = size
    
    def store_dictio(self, dictio, name_dictio=None, script=None):
        if isinstance(dictio, dict) and name_dictio is None: 
            raise Exception(f'The first parameter is not a dictionary, you must write a name for it (second parameter).')
        if dictio is not None and dictio: # Store if it is not None or empty
            if script is None: 

                dictionary = self.__stats
            else:                
                if script not in self.__stats:
                    self.__stats[script] = {} 
                dictionary = self.__stats[script]

            if name_dictio is not None:
                if name_dictio in dictionary:
                    raise Exception(f'The key \"{name_data}\" contains data, choose another name or remove the key first.')
                else:
                    dictionary[name_dictio] = dictio
            elif not any(key in self.__stats for key in dictio):
                dictionary.update(dictio)   
    
    def delete_stat(self, stat):
        if stat in self.__stats:
            stat_data = deepcopy(self.__stats.get(stat))
            del self.__stats[stat]
            return stat_data
        else:
            print("WARNING: The state that was attempted to be deleted does not exist.")
            return None

    def get_stats(self, stats=None):
        if stats is None:
            return deepcopy(self.__stats)
        
        dictio = deepcopy(self.__stats.get(stats))
        if dictio is None:
            raise Exception(f'The key \"{stats}\" not exist.')
        return dictio

    def load_stats_to_json(self, bucket, name_file='log_stats'):
        storage = Storage()
        storage.put_object(bucket=bucket, key=f'stats/{name_file}.json', body=str(json.dumps(self.__stats,indent=2)))    