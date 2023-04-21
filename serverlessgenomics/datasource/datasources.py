from enum import Enum, auto


class FASTQSource(Enum):
    S3_GZIP = auto()
    SRA = auto()


class FASTASource(Enum):
    S3_FASTA = auto()
