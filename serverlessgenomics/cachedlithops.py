from __future__ import annotations

import logging
import os.path
from typing import TYPE_CHECKING, Optional, List, Union, Tuple, Dict, Any, Callable

import lithops
import diskcache
from lithops.future import ResponseFuture

if TYPE_CHECKING:
    from .parameters import PipelineRun

CACHE_PATH = os.path.expanduser(os.path.join('~', 'serverless-genomics.cache'))

logger = logging.getLogger(__name__)


class CachedLithopsInvoker:
    def __init__(self, pipeline_params: PipelineRun, **kwargs):
        self.__fexec = lithops.FunctionExecutor(**kwargs)
        self.__cache = diskcache.Cache(CACHE_PATH)
        self.pipeline_params = pipeline_params

    def _get_cache_key(self, key):
        return '_'.join([self.pipeline_params.execution_id, key])

    def clear_all_cache(self):
        self.__cache.clear()

    def clear_cache_run(self):
        keys = filter(lambda key: key.strartswith(self.pipeline_params.execution_id), self.__cache.iterkeys())
        for key in keys:
            self.__cache.delete(key)

    def call_async(self, func, data, extra_env=None, runtime_memory=None, timeout=None, include_modules=[],
                   exclude_modules=[]):
        key = self._get_cache_key(func.__name__)
        if key in self.__cache:
            logger.debug('Getting futures from cache (key=%s)', key)
            fut = self.__cache[key]
        else:
            fut = self.__fexec.call_async(func, data, extra_env, runtime_memory, timeout,
                                          include_modules, exclude_modules)
            self.__cache[self._get_cache_key(func.__name__)] = fut
            logger.debug('Stored futures to cache (key=%s)', key)
        return fut

    def call(self, func, data, extra_env=None, runtime_memory=None, timeout=None, include_modules=[],
             exclude_modules=[]):
        key = self._get_cache_key(func.__name__)
        if key in self.__cache:
            logger.debug('Getting futures from cache (key=%s)', key)
            fut = self.__cache[key]
        else:
            fut = self.__fexec.call_async(func, data, extra_env, runtime_memory, timeout,
                                          include_modules, exclude_modules)
            self.__cache[key] = fut
            logger.debug('Stored futures to cache (key=%s)', key)
        result = self.__fexec.get_result(fs=fut)
        return result

    def map_async(self, map_function, map_iterdata, chunksize=None, extra_args=None, extra_env=None,
                  runtime_memory=None, obj_chunk_size=None, obj_chunk_number=None, obj_newline='\n', timeout=None,
                  include_modules=[], exclude_modules=[]):
        raise NotImplementedError()

    def map(self, map_function, map_iterdata, chunksize=None, extra_args=None, extra_env=None, runtime_memory=None,
            obj_chunk_size=None, obj_chunk_number=None, obj_newline='\n', timeout=None, include_modules=[],
            exclude_modules=[]):
        key = self._get_cache_key(map_function.__name__)
        if key in self.__cache:
            logger.debug('Getting futures from cache (key=%s)', key)
            fut = self.__cache[key]
        else:
            fut = self.__fexec.map(map_function, map_iterdata, chunksize, extra_args, extra_env, runtime_memory,
                                   obj_chunk_size, obj_chunk_number, obj_newline, timeout, include_modules,
                                   exclude_modules)
            self.__cache[key] = fut
            logger.debug('Stored futures to cache (key=%s)', key)
        result = self.__fexec.get_result(fs=fut)
        return result

    def map_reduce_async(self, map_function, map_iterdata, reduce_function, chunksize=None, extra_args=None,
                         extra_env=None, map_runtime_memory=None, reduce_runtime_memory=None, timeout=None,
                         obj_chunk_size=None, obj_chunk_number=None, obj_newline='\n', obj_reduce_by_key=None,
                         spawn_reducer=20, include_modules=[], exclude_modules=[]):
        raise NotImplementedError()

    def map_reduce(self, map_function, map_iterdata, reduce_function, chunksize=None, extra_args=None, extra_env=None,
                   map_runtime_memory=None, reduce_runtime_memory=None, timeout=None, obj_chunk_size=None,
                   obj_chunk_number=None, obj_newline='\n', obj_reduce_by_key=None, spawn_reducer=20,
                   include_modules=[], exclude_modules=[]):
        key = self._get_cache_key(map_function.__name__ + '-' + reduce_function.__name__)
        if key in self.__cache:
            logger.debug('Getting futures from cache (key=%s)', key)
            fut = self.__cache[key]
        else:
            fut = self.__fexec.map_reduce(map_function, map_iterdata, reduce_function, chunksize, extra_args, extra_env,
                                          map_runtime_memory, reduce_runtime_memory, timeout, obj_chunk_size,
                                          obj_chunk_number, obj_newline, obj_reduce_by_key, spawn_reducer,
                                          include_modules, exclude_modules)
            self.__cache[key] = fut
            logger.debug('Stored futures to cache (key=%s)', key)
        result = self.__fexec.get_result(fs=fut)
        return result
