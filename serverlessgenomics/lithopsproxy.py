from __future__ import annotations
from typing import TYPE_CHECKING, Optional, List, Union, Tuple, Dict, Any, Callable

import lithops
from lithops.future import ResponseFuture

if TYPE_CHECKING:
    from .parameters import PipelineRun


class LithopsProxy:
    def __init__(self,
                 mode: Optional[str] = None,
                 config: Optional[Dict[str, Any]] = None,
                 backend: Optional[str] = None,
                 storage: Optional[str] = None,
                 runtime: Optional[str] = None,
                 runtime_memory: Optional[int] = None,
                 monitoring: Optional[str] = None,
                 max_workers: Optional[int] = None,
                 worker_processes: Optional[int] = None,
                 remote_invoker: Optional[bool] = None,
                 log_level: Optional[str] = False):
        self.__fexec = lithops.FunctionExecutor(mode, config, backend, storage, runtime, runtime_memory, monitoring,
                                                max_workers, worker_processes, remote_invoker, log_level)
        self.__futures = {}

    def call_async(self,
                   stage: str,
                   pipeline_params: PipelineRun,
                   func: Callable,
                   data: Union[List[Any], Tuple[Any, ...], Dict[str, Any]],
                   extra_env: Optional[Dict] = None,
                   runtime_memory: Optional[int] = None,
                   timeout: Optional[int] = None,
                   include_modules: Optional[List] = [],
                   exclude_modules: Optional[List] = []) -> ResponseFuture:
        fut = self.__fexec.call_async(func, data, extra_env, runtime_memory,
                                      timeout, include_modules, exclude_modules)
        self.__futures[stage] = fut
        return fut

    def call(self,
             stage: str,
             pipeline_params: PipelineRun,
             func: Callable,
             data: Union[List[Any], Tuple[Any, ...], Dict[str, Any]],
             extra_env: Optional[Dict] = None,
             runtime_memory: Optional[int] = None,
             timeout: Optional[int] = None,
             include_modules: Optional[List] = [],
             exclude_modules: Optional[List] = []):
        fut = self.__fexec.call_async(func, data, extra_env, runtime_memory,
                                      timeout, include_modules, exclude_modules)
        self.__futures[stage] = fut
        result = self.__fexec.get_result(fs=fut)
        return result
