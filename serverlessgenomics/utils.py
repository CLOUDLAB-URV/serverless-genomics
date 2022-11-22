import logging
import json
from .parameters import PipelineParameters
from dataclasses import asdict

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
