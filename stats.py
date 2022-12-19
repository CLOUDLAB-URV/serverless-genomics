'''
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
     __stats.timer('script1').start()
     ...
     __stats.timer('script1').stop()

     __stats.value('data_size', 100)
     return ..., __stats


TUPLES?
'''
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
        else:
            self.__tmp_registrer[script] = {} 

        self.__tmp_registrer[script]["start_time"] = time.time()
    
    def timer_stop(self, script):
        end_time = time.time()
        
        if script in self.__tmp_registrer and "start_time" in self.__tmp_registrer[script]:
            if script not in self.__stats:
                self.__stats[script] = {} 
            self.__stats[script]["execution_time"] = end_time - self.__tmp_registrer[script]["start_time"]
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

    def load_stats_to_json(self, name_file="logs_stats"):
        with open(f'{name_file}.txt', 'w') as json_file:
            json.dump(self.__stats, json_file)



    