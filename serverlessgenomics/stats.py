'''
TODO
    1. Crear parametro para decir si se quiere almacenar los stats o no
    2. Poner timers en las funciones que necesarias
    3. Recoger tamaños de ficheros que se traspasan
        - /Index fastq
        - Index fasta
        - Fasta
per exemple, pots tindre una classe Stats, on puguis registrar timers i valors. 
Durant la funció lambda vas registrant timers (temps d'execució de scripts i les coses de la funcion) 
i valors (mides de fixters principalment) i despres retornar aquesta classe juntament amb els resultats


despres al client si s'indica per parametre, ho pot guardar en json o alguna cosa aixi

si fas un lithops.map de my_function, tindras una llista dels resultats que son tuples amb el __stats
10:33
s'hauria de treure de les tuples per nomes tindre una llista de __stats, i tindre despres, per exemple, 
una funcio stats_to_json(__stats) que rep una llista de __stats i que ho passa a json


sobre les funcions que es criden amb lithops principalment
__tmp_registrer

ef my_function():
     __stats = Stats()
     __stats.perf_counterr('script1').start()
     ...
     __stats.perf_counterr('script1').stop()

     __stats.value('data_size', 100)
     return ..., __stats


TUPLES?

penso que estaria be tindre uns stats per cada stage del pipeline, 
indicant el temps d'execució de cada script, mida de les dades intermitjes, 
quan triga en pujar / descarregar de s3... tot aixo
'''
import storage
import time
import json
from copy import deepcopy

class Stats:
    def __init__(self):
        self.__tmp_registrer={}
        self.__stats={}

    def timer_start(self, script):
        if script in self.__tmp_registrer:
            print(f'WARNING: the counter of the timer \"{script}\" was already running, it will restart.')
        elif script in self.__stats and "execution_time" in self.__stats[script]:
            raise Exception(f'The timer of \"{script}\" already existed, choose another name or remove the existing timer first.')
 
        self.__tmp_registrer[script] = time.perf_counter()
    
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
        if name_data in self.__stats and script == None:
            raise Exception(f'The key \"{name_data}\" contains data, choose another name or remove the key first.')
        elif script is None: 
            self.__stats[name_data] = size
        else:
            if script not in self.__stats:
                self.__stats[script] = {} 
            self.__stats[script][name_data] = size
    
    def store_dictio(self, dictio, name_dictio=None):
        if not isinstance(dictio, dict): 
            raise Exception(f'The second parameter must be a dictionary.')

        if name_dictio is not None:
            if name_dictio in self.__stats:
                raise Exception(f'The key \"{name_data}\" contains data, choose another name or remove the key first.')
            else:
               self.__stats[name_dictio] = dictio
        elif not any(key in self.__stats for key in dictio):
            self.__stats.update(dictio)

        


            
    
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

    def load_stats_to_json(self, bucket, name_file=):
        put_object(bucket=bucket, key=f'stats/{name_file}.json', body=str(delf.__stats))    