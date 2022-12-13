from __future__ import annotations

import logging
import os.path
from typing import TYPE_CHECKING

import lithops
import shelve
import os

if TYPE_CHECKING:
    from .parameters import PipelineRun

from .constants import CACHE_PATH


logger = logging.getLogger(__name__)


class CachedLithopsInvoker:
    def __init__(self, pipeline_params: PipelineRun, **kwargs):
        self.__fexec = lithops.FunctionExecutor(**kwargs)
        self.pipeline_params = pipeline_params

    def _fmt_full_key(self, key):
        return os.path.join(self.pipeline_params.run_id, 'lithops', key)

    def _clear_all_cache(self):
        with shelve.open(CACHE_PATH) as cache:
            cache.clear()

    def _clear_cache_run(self):
        with shelve.open(CACHE_PATH) as cache:
            keys = filter(lambda key: key.strartswith(self.pipeline_params.run_id), cache.keys())
            for key in keys:
                cache.pop(key)

    def _get_cache(self, key):
        full_key = self._fmt_full_key(key)
        with shelve.open(CACHE_PATH) as cache:
            if full_key in cache:
                logger.debug('Got futures from cache (key=%s)', full_key)
                return cache[full_key]
            else:
                return None

    def _save_cache(self, key, futures):
        full_key = self._fmt_full_key(key)
        with shelve.open(CACHE_PATH) as cache:
            cache[full_key] = futures
        logger.debug('Stored futures to cache (key=%s)', full_key)

    def call_async(self, func, data, extra_env=None, runtime_memory=None, timeout=None, include_modules=[],
                   exclude_modules=[]):
        fut = self._get_cache(func.__name__)
        if fut is None:
            fut = self.__fexec.call_async(func, data, extra_env, runtime_memory, timeout,
                                          include_modules, exclude_modules)
            # self._save_cache(func.__name__, fut)
        return fut

    def call(self, func, data, extra_env=None, runtime_memory=None, timeout=None, include_modules=[],
             exclude_modules=[]):
        fut = self._get_cache(func.__name__)
        if fut is None:
            fut = self.__fexec.call_async(func, data, extra_env, runtime_memory, timeout,
                                          include_modules, exclude_modules)
            # self._save_cache(func.__name__, fut)
        result = self.__fexec.get_result(fs=fut)
        self._save_cache(func.__name__, fut)
        return result

    def map_async(self, map_function, map_iterdata, chunksize=None, extra_args=None, extra_env=None,
                  runtime_memory=None, obj_chunk_size=None, obj_chunk_number=None, obj_newline='\n', timeout=None,
                  include_modules=[], exclude_modules=[]):
        fut = self._get_cache(map_function.__name__)
        if fut is None:
            fut = self.__fexec.map(map_function, map_iterdata, chunksize, extra_args, extra_env, runtime_memory,
                                   obj_chunk_size, obj_chunk_number, obj_newline, timeout, include_modules,
                                   exclude_modules)
            # self._save_cache(map_function.__name__, fut)
        return fut

    def map(self, map_function, map_iterdata, chunksize=None, extra_args=None, extra_env=None, runtime_memory=None,
            obj_chunk_size=None, obj_chunk_number=None, obj_newline='\n', timeout=None, include_modules=[],
            exclude_modules=[]):
        fut = self._get_cache(map_function.__name__)
        if fut is None:
            fut = self.__fexec.map(map_function, map_iterdata, chunksize, extra_args, extra_env, runtime_memory,
                                   obj_chunk_size, obj_chunk_number, obj_newline, timeout, include_modules,
                                   exclude_modules)
            # self._save_cache(map_function.__name__, fut)
        result = self.__fexec.get_result(fs=fut)
        self._save_cache(map_function.__name__, fut)
        return result

    def map_reduce_async(self, map_function, map_iterdata, reduce_function, chunksize=None, extra_args=None,
                         extra_env=None, map_runtime_memory=None, reduce_runtime_memory=None, timeout=None,
                         obj_chunk_size=None, obj_chunk_number=None, obj_newline='\n', obj_reduce_by_key=None,
                         spawn_reducer=20, include_modules=[], exclude_modules=[]):
        fut = self._get_cache(map_function.__name__ + '-' + reduce_function.__name__)
        if fut is None:
            fut = self.__fexec.map_reduce(map_function, map_iterdata, reduce_function, chunksize, extra_args, extra_env,
                                          map_runtime_memory, reduce_runtime_memory, timeout, obj_chunk_size,
                                          obj_chunk_number, obj_newline, obj_reduce_by_key, spawn_reducer,
                                          include_modules, exclude_modules)
            # self._save_cache(map_function.__name__ + '-' + reduce_function.__name__, fut)
        return fut

    def map_reduce(self, map_function, map_iterdata, reduce_function, chunksize=None, extra_args=None, extra_env=None,
                   map_runtime_memory=None, reduce_runtime_memory=None, timeout=None, obj_chunk_size=None,
                   obj_chunk_number=None, obj_newline='\n', obj_reduce_by_key=None, spawn_reducer=20,
                   include_modules=[], exclude_modules=[]):
        fut = self._get_cache(map_function.__name__ + '-' + reduce_function.__name__)
        if fut is None:
            fut = self.__fexec.map_reduce(map_function, map_iterdata, reduce_function, chunksize, extra_args, extra_env,
                                          map_runtime_memory, reduce_runtime_memory, timeout, obj_chunk_size,
                                          obj_chunk_number, obj_newline, obj_reduce_by_key, spawn_reducer,
                                          include_modules, exclude_modules)
            # self._save_cache(map_function.__name__ + '-' + reduce_function.__name__, fut)
        result = self.__fexec.get_result(fs=fut)
        self._save_cache(map_function.__name__ + '-' + reduce_function.__name__, fut)
        return result
