import logging
import re
import subprocess
import tempfile
import time
import threading
import os
import io


import lithops

from ..utils import get_gztool_path, S3Path, force_delete_local_path

logger = logging.getLogger(__name__)

CHUNK_SIZE = 65536


def fetch_fastq_chunk(
    fastq_chunk: dict,
    target_filename: str,
    storage: lithops.Storage,
    fastq_path: S3Path,
    storage_bucket: str,
    fastqgz_index_key: str,
):
    tmp_index_file = tempfile.mktemp()
    gztool = get_gztool_path()
    lines = []
    lines_to_read = fastq_chunk["line_1"] - fastq_chunk["line_0"] + 1

    try:
        t0 = time.perf_counter()

        # Get index and store it to temp file
        storage.download_file(bucket=storage_bucket, key=fastqgz_index_key, file_name=tmp_index_file)

        # Get compressed byte range
        extra_get_args = {"Range": f"bytes={fastq_chunk['range_0'] - 1}-{fastq_chunk['range_1'] - 1}"}
        body = storage.get_object(fastq_path.bucket, fastq_path.key, True, extra_get_args)

        cmd = [gztool, "-I", tmp_index_file, "-n", str(fastq_chunk["range_0"]), "-L", str(fastq_chunk["line_0"])]
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

        # TODO program might get stuck if subprocess fails, blocking io should be done in a backgroun thread or using async/await
        def _writer_feeder():
            logger.debug("Writer thread started")
            input_chunk = body.read(CHUNK_SIZE)
            while input_chunk != b"":
                # logger.debug('Writing %d bytes to pipe', len(chunk))
                try:
                    proc.stdin.write(input_chunk)
                except BrokenPipeError:
                    break
                input_chunk = body.read(CHUNK_SIZE)
            try:
                proc.stdin.flush()
                proc.stdin.close()
            except BrokenPipeError:
                pass
            logger.debug("Writer thread finished")

        writer_thread = threading.Thread(target=_writer_feeder)
        writer_thread.start()

        output_chunk = proc.stdout.read(CHUNK_SIZE)
        last_line = None
        while output_chunk != b"":
            # logger.debug('Read %d bytes from pipe', len(chunk))
            text = output_chunk.decode("utf-8")
            chunk_lines = text.splitlines()
            if last_line is not None:
                last_line = last_line + chunk_lines.pop(0)
                lines.append(last_line)
                last_line = None
            if text[-1] != "\n":
                last_line = chunk_lines.pop()

            lines.extend(chunk_lines)

            # Stop decompressing lines if number of lines to read in this chunk is reached
            if len(lines) > lines_to_read:
                proc.stdout.close()
                break

            # Try to read next decompressed chunk
            # a ValueError is raised if the pipe is closed, meaning the writer or the subprocess closed it
            try:
                output_chunk = proc.stdout.read(CHUNK_SIZE)
            except ValueError:
                output_chunk = b""

        try:
            proc.wait()
        except ValueError as e:
            logger.error(e)

        writer_thread.join()

        t1 = time.perf_counter()
        logger.debug("Got partition in %.3f seconds", t1 - t0)
        # TODO write lines to file as decompressed instead of saving them all in memory
        with open(target_filename, "w") as target_file:
            target_file.writelines((line + "\n" for line in lines[: fastq_chunk["line_1"] - fastq_chunk["line_0"]]))

    finally:
        force_delete_local_path(tmp_index_file)


def fetch_fasta_chunk(fasta_chunk: dict, target_filename: str, storage: lithops.Storage, fasta_path: S3Path):
    # Get header data
    extra_args = {"Range": f"bytes={fasta_chunk['offset_head']}-{fasta_chunk['offset_base']}"}
    header_body = storage.get_object(bucket=fasta_path.bucket, key=fasta_path.key, extra_get_args=extra_args)
    chunk_body = list(re.finditer(r">.+\n", header_body.decode("utf-8")))[0].group()

    # Get chunk body and append to header
    extra_args = {"Range": f"bytes={fasta_chunk['offset_base']}-{fasta_chunk['last_byte']}"}
    base = storage.get_object(bucket=fasta_path.bucket, key=fasta_path.key, extra_get_args=extra_args).decode("utf-8")
    chunk_body += base[1::] if base[0:1] == "\n" else base

    with open(target_filename, "w") as target_file:
        target_file.writelines(chunk_body)


def fetch_fastq_chunk_sra(seq_name: str, fastq_chunk: dict, target_filename: str):
    """
    Function to retrieve the relevant SRA chunk using fastq-dump and save it to object storage
    """

    start_read = int(fastq_chunk["read_0"])
    end_read = int(fastq_chunk["read_1"])

    # To suppress a warning that appears the first time vdb-config is used
    subprocess.run(["vdb-config", "-i"])
    # Report cloud identity so it can take data from s3 needed to be executed only once per vm
    subprocess.run(["vdb-config", "--report-cloud-identity", "yes"], capture_output=True)

    # Run fastq-dump with the specified range of reads, splits files in two files if paired end
    subprocess.run(
        ["fastq-dump", "--split-files", seq_name, "-X", str(start_read), "-N", str(end_read)], capture_output=True
    )

    original_output_file_1 = f"{seq_name}_1.fastq"
    new_output_file_1 = f"{target_filename}"

    os.rename(original_output_file_1, new_output_file_1)

    # Save the output file to the desired target filename in object storage

    print(f"Finished fetching chunk {fastq_chunk['chunk_id']} and saved to object storage")
