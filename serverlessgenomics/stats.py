from __future__ import annotations

import logging

from lithops import Storage
import time
import json
from copy import deepcopy

logger = logging.getLogger(__name__)


class _ContextManagerTimer:
    def __init__(self, key: str, stats: Stats):
        self.__key = key
        self.__stats = stats

    def __enter__(self):
        self.__stats.start_timer(self.__key)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__stats.stop_timer(self.__key)


class Stats:
    def __init__(self):
        self.__tmp_registrer = {}
        self.__stats = {}

        self.__timers = {}
        self.__values = {}

    def start_timer(self, key):
        if key in self.__timers:
            logger.warning("Timer %s was already running, it will restart")
        self.__timers[key] = {'t0': time.time(), 't0_perf_counter': time.perf_counter()}

    def stop_timer(self, key):
        if key not in self.__timers:
            logger.warning("Timer %s not registered, skipping...")
            return

        self.__timers[key]['t1'] = time.time()
        self.__timers[key]['elapsed'] = time.perf_counter() - self.__timers[key]['t0_perf_counter']
        del self.__timers[key]['t0_perf_counter']

    def set_value(self, key, value):
        if key in self.__values:
            logger.warning("Value with key %s already exists, it will overwrite it")
        self.__values[key] = value

    def incr_value(self, key, delta=1):
        if key not in self.__values:
            self.__values[key] = delta
        else:
            self.__values[key] += delta

    def timeit(self, key):
        return _ContextManagerTimer(key, self)

    def timer_start(self, script, extra_time=None):
        if script in self.__tmp_registrer:
            print(f'WARNING: the counter of the timer \"{script}\" was already running, it will restart.')
        elif script in self.__stats and "execution_time" in self.__stats[script]:
            raise Exception(
                f'The timer of \"{script}\" already existed, choose another name or remove the existing timer first.')

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
            raise Exception(
                f'The first parameter is not a dictionary, you must write a name for it (second parameter).')
        if dictio is not None and dictio:  # Store if it is not None or empty
            if script is None:

                dictionary = self.__stats
            else:
                if script not in self.__stats:
                    self.__stats[script] = {}
                dictionary = self.__stats[script]

            if name_dictio is not None:
                if name_dictio in dictionary:
                    raise Exception(
                        f'The key \"{name_dictio}\" contains data, choose another name or remove the key first.')
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
        storage.put_object(bucket=bucket, key=f'stats/{name_file}.json', body=str(json.dumps(self.__stats, indent=2)))


def lithops_function(f):
    def _wrapper():
        print(f'Starting execution of {f.__name__}')
        stats = Stats()
        with stats.timeit('function'):
            res = f()
        print('Execution of {} ended - took {:.2f}'.format(f.__name__, stats.timers['function']))
        assert isinstance(res, dict)
        res['stats'] = stats
        return res
    return _wrapper
