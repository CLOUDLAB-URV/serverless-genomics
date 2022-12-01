import logging
import json
from .parameters import PipelineParameters
from dataclasses import asdict
from .preprocessing import create_fasta_chunk_for_runtime
from lithops import Storage
import shutil

logger = logging.getLogger(__name__)


def setup_logging(level=logging.INFO):
    root_logger = logging.getLogger('serverlessgenomics')
    root_logger.setLevel(level)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s")
    ch.setFormatter(formatter)
    root_logger.addHandler(ch)


def log_parameters(params: PipelineParameters):
    logger.debug('Pipeline parameters: \n %s', json.dumps(asdict(params), indent=2))


def copy_to_runtime(storage: Storage, bucket: str, folder: str, file_name: str, byte_range={}, fasta=None):
    """
    Copy file from S3 to local storage.
    """
    temp_file = "/tmp/" + file_name
    with open(temp_file, 'wb') as file:
        if fasta is not None: 
            file.write(create_fasta_chunk_for_runtime(storage, bucket, fasta, byte_range, folder, file_name))
        else:
            shutil.copyfileobj(storage.get_object(bucket=bucket, key=folder+file_name, stream=True, extra_get_args=byte_range), file)
    return temp_file
