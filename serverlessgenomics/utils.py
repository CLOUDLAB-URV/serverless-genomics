from __future__ import annotations

import logging
import json
import os
import shutil
import subprocess

from contextlib import suppress
from pathlib import PurePath, _PosixFlavour
from typing import TYPE_CHECKING
from dataclasses import asdict

import lithops
from lithops.storage.utils import StorageNoSuchKeyError

# from .preprocessing.preprocess_fasta import create_fasta_chunk_for_runtime

if TYPE_CHECKING:
    from lithops import Storage
    from .parameters import PipelineRun

logger = logging.getLogger(__name__)


class _S3Flavour(_PosixFlavour):
    is_supported = True

    def parse_parts(self, parts):
        drv, root, parsed = super().parse_parts(parts)
        for part in parsed[1:]:
            if part == "..":
                index = parsed.index(part)
                parsed.pop(index - 1)
                parsed.remove(part)
        return drv, root, parsed

    def make_uri(self, path):
        uri = super().make_uri(path)
        return uri.replace("file:///", "s3://")


class S3Path(PurePath):
    """
    PurePath subclass for AWS S3 service.
    Source: https://github.com/liormizr/s3path
    S3 is not a file-system, but we can look at it like a POSIX system.
    """

    _flavour = _S3Flavour()
    __slots__ = ()

    @classmethod
    def from_uri(cls, uri: str) -> "S3Path":
        """
        from_uri class method create a class instance from url

        >> from s3path import PureS3Path
        >> PureS3Path.from_url('s3://<bucket>/<key>')
        << PureS3Path('/<bucket>/<key>')
        """
        if not uri.startswith("s3://"):
            raise ValueError("Provided uri seems to be no S3 URI!")
        return cls(uri[4:])

    @classmethod
    def from_bucket_key(cls, bucket: str, key: str) -> "S3Path":
        """
        from_bucket_key class method create a class instance from bucket, key pair's

        >> from s3path import PureS3Path
        >> PureS3Path.from_bucket_key(bucket='<bucket>', key='<key>')
        << PureS3Path('/<bucket>/<key>')
        """
        bucket = cls(cls._flavour.sep, bucket)
        if len(bucket.parts) != 2:
            raise ValueError("bucket argument contains more then one path element: {}".format(bucket))
        key = cls(key)
        if key.is_absolute():
            key = key.relative_to("/")
        return bucket / key

    @property
    def bucket(self) -> str:
        """
        The AWS S3 Bucket name, or ''
        """
        self._absolute_path_validation()
        with suppress(ValueError):
            _, bucket, *_ = self.parts
            return bucket
        return ""

    @property
    def key(self) -> str:
        """
        The AWS S3 Key name, or ''
        """
        self._absolute_path_validation()
        key = self._flavour.sep.join(self.parts[2:])
        return key

    @property
    def virtual_directory(self) -> str:
        """
        The parent virtual directory of a key
        Example: foo/bar/baz -> foo/baz
        """
        vdir, _ = self.key.rsplit("/", 1)
        return vdir

    def as_uri(self) -> str:
        """
        Return the path as a 's3' URI.
        """
        return super().as_uri()

    def _absolute_path_validation(self):
        if not self.is_absolute():
            raise ValueError("relative path have no bucket, key specification")

    def __repr__(self) -> str:
        return "{}(bucket={},key={})".format(self.__class__.__name__, self.bucket, self.key)


def force_delete_local_path(path):
    if os.path.exists(path):
        if os.path.isfile(path):
            os.remove(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)


def try_head_object(storage: lithops.Storage, bucket: str, key: str):
    try:
        res = storage.head_object(bucket, key)
        return res
    except StorageNoSuchKeyError:
        return None


def try_get_object(storage: lithops.Storage, bucket: str, key: str, stream: bool = False, extra_get_args: dict = None):
    try:
        res = storage.get_object(bucket, key, stream, extra_get_args)
        return res
    except StorageNoSuchKeyError:
        return None


def setup_logging(level=logging.INFO):
    root_logger = logging.getLogger("serverlessgenomics")
    root_logger.setLevel(level)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(name)s.%(funcName)s:%(lineno)d -- %(message)s")
    ch.setFormatter(formatter)
    root_logger.addHandler(ch)


def log_parameters(params: PipelineRun):
    for k, v in asdict(params).items():
        logger.debug("\t\t%s = %s", k, repr(v))


def get_gztool_path():
    """
    Utility function that returns the absolute path for gzip file binary or raises exception if it is not found
    """
    proc = subprocess.run(["which", "gztool"], check=True, capture_output=True, text=True)
    path = proc.stdout.rstrip("\n")
    logger.debug("Using gztool located in %s", path)
    return path


def copy_to_runtime(storage: Storage, bucket: str, folder: str, file_name: str, byte_range={}, fasta=None):
    """
    Copy file from S3 to local storage.
    """
    temp_file = "/tmp/" + file_name
    with open(temp_file, "wb") as file:
        if fasta is not None:
            # file.write(create_fasta_chunk_for_runtime(storage, bucket, fasta, byte_range, folder, file_name))
            raise Exception()
        else:
            shutil.copyfileobj(
                storage.get_object(bucket=bucket, key=folder + file_name, stream=True, extra_get_args=byte_range), file
            )
    return temp_file


def split_data_result(result):
    aux_timer = []
    aux_info = []
    for info, timer in result:
        aux_info.append(info)
        aux_timer.append(timer)
    return tuple(aux_info), aux_timer
