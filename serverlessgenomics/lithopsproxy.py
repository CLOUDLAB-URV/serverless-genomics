from __future__ import annotations
from typing import TYPE_CHECKING, Optional, List, Union, Tuple, Dict, Any, Callable

import lithops
from lithops.future import ResponseFuture

if TYPE_CHECKING:
    from .parameters import PipelineRun


class LithopsProxy:
    def __init__(self, **kwargs):
        self.__fexec = lithops.FunctionExecutor(**kwargs)
        self.__futures = {}

    def call_async(self, stage: str, pipeline_params: PipelineRun, *args, **kwargs) -> ResponseFuture:
        fut = self.__fexec.call_async(*args, **kwargs)
        self.__futures[stage] = fut
        return fut

    def call(self, stage: str, pipeline_params: PipelineRun, *args, **kwargs):
        fut = self.__fexec.call_async(*args, **kwargs)
        self.__futures[stage] = fut
        result = self.__fexec.get_result(fs=fut)
        return result

    def map_async(self, stage, pipeline_params, *args, **kwargs):
        raise NotImplementedError()

    def map(self, stage, pipeline_params, *args, **kwargs):
        raise NotImplementedError()

    def map_reduce_async(self, stage, pipeline_params, *args, **kwargs):
        raise NotImplementedError()

    def map_reduce(self, stage, pipeline_params, *args, **kwargs):
        fut = self.__fexec.map_reduce(*args, **kwargs)
        self.__futures[stage] = fut
        result = self.__fexec.get_result(fs=fut)
        return result
