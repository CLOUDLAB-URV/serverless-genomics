from __future__ import annotations

import logging

import lithops

logger = logging.getLogger(__name__)


class LithopsInvokerWrapper:
    def __init__(self, lithops_config: dict):
        config = lithops_config or {}
        self.__fexec = lithops.FunctionExecutor(**config)

    def call(
        self, func, data, /, extra_env=None, runtime_memory=None, timeout=None, include_modules=[], exclude_modules=[]
    ):
        fut = self.__fexec.call_async(func, data, extra_env, runtime_memory, timeout, include_modules, exclude_modules)
        result = self.__fexec.get_result(fs=fut)
        return result

    def map(
        self,
        map_function,
        map_iterdata,
        chunksize=None,
        extra_args=None,
        extra_env=None,
        runtime_memory=None,
        obj_chunk_size=None,
        obj_chunk_number=None,
        obj_newline="\n",
        timeout=None,
        include_modules=[],
        exclude_modules=[],
    ):
        fut = self.__fexec.map(
            map_function,
            map_iterdata,
            chunksize,
            extra_args,
            extra_env,
            runtime_memory,
            obj_chunk_size,
            obj_chunk_number,
            obj_newline,
            timeout,
            include_modules,
            exclude_modules,
        )
        res = self.__fexec.get_result(fs=fut)
        return res

    def map_reduce(
        self,
        map_function,
        map_iterdata,
        reduce_function,
        chunksize=None,
        extra_args=None,
        extra_args_reduce=None,
        extra_env=None,
        map_runtime_memory=None,
        reduce_runtime_memory=None,
        timeout=None,
        obj_chunk_size=None,
        obj_chunk_number=None,
        obj_newline="\n",
        obj_reduce_by_key=False,
        spawn_reducer=20,
        include_modules=[],
        exclude_modules=[],
    ):
        fut = self.__fexec.map_reduce(
            map_function,
            map_iterdata,
            reduce_function,
            chunksize,
            extra_args,
            extra_args_reduce,
            extra_env,
            map_runtime_memory,
            reduce_runtime_memory,
            timeout,
            obj_chunk_size,
            obj_chunk_number,
            obj_newline,
            obj_reduce_by_key,
            spawn_reducer,
            include_modules,
            exclude_modules,
        )
        res = self.__fexec.get_result(fs=fut)
        return res
