import os
import tempfile
import argparse
import io

import asyncio
from aiobotocore import session as aiobotocore_sess

import lithops


STORAGE_PROGRAMS_FOLDER = 'programs/'

# https://github.com/aitorarjona/active-s3-modules/blob/f4e09edd0254d5b2dc0470431d802fd67019baf6/lidar/app.py#L27


class MultipartUploader:
    PART_SIZE = 64 * 1024 * 1024  # 64 MiB

    def __init__(self, s3client, bucket, key, content_type="application/octet-stream"):
        self.__buff = io.BytesIO()
        self.__part_count = 1
        self.__s3client = s3client
        self.__bucket = bucket
        self.__key = key
        self.__closed = False
        self.__parts = []
        self._content_type = content_type
        self.etag = None

    async def setup(self):
        multipart_res = await self.__s3client.create_multipart_upload(
            Bucket=self.__bucket,
            Key=self.__key,
            ContentType=self._content_type
        )
        print(multipart_res)

        self.__upload_id = multipart_res['UploadId']

    async def write(self, data):
        self.__buff.write(data)
        buff_sz = self.__buff.tell()
        #print(f'buff size is {buff_sz}')
        if buff_sz >= self.PART_SIZE:
            self.__buff.seek(0)
            chunk = self.__buff.read(self.PART_SIZE)
            upload_part_res = await self.__s3client.upload_part(
                Body=chunk,
                Bucket=self.__bucket,
                Key=self.__key,
                PartNumber=self.__part_count,
                UploadId=self.__upload_id
            )
            print(upload_part_res)
            self.__parts.append({
                'PartNumber': self.__part_count,
                'ETag': upload_part_res['ETag']
            })
            self.__part_count += 1
            rest = self.__buff.read(buff_sz - self.PART_SIZE)
            self.__buff = io.BytesIO()
            self.__buff.write(rest)

    async def flush(self):
        print('flushing buffer')
        if self.__buff.tell() > 0:
            self.__buff.seek(0)
            chunk = self.__buff.read(self.PART_SIZE)
            print(f'put part of size {len(chunk)}')
            while len(chunk) > 0:
                upload_part_res = await self.__s3client.upload_part(
                    Body=chunk,
                    Bucket=self.__bucket,
                    Key=self.__key,
                    PartNumber=self.__part_count,
                    UploadId=self.__upload_id
                )
                print(upload_part_res)
                self.__parts.append({
                    'PartNumber': self.__part_count,
                    'ETag': upload_part_res['ETag']
                })
                self.__part_count += 1
                chunk = self.__buff.read(self.PART_SIZE)

        print(self.__parts)

        complete_multipart_response = await self.__s3client.complete_multipart_upload(
            Bucket=self.__bucket,
            Key=self.__key,
            UploadId=self.__upload_id,
            MultipartUpload={'Parts': self.__parts}
        )
        print(complete_multipart_response)

        self.etag = complete_multipart_response['ETag']

    async def abort_upload(self):
        if self.__upload_id:
            print('aborting upload')
            abort_response = await self.__s3client.abort_multipart_upload(
                Bucket=self.__bucket,
                Key=self.__key,
                UploadId=self.__upload_id
            )
            print(abort_response)


def download_binary(storage, args):
    """
    Alternatively create a custom docker runtime in a private docker
    registry and save the binary to any folder
    """
    binary_name = args.program_and_args[0][2:]  # discard "./"

    tmp_prefixed_binary_name = f"/tmp/{binary_name}"

    storage.download_file(storage.bucket,
                          STORAGE_PROGRAMS_FOLDER + binary_name, tmp_prefixed_binary_name)

    os.chmod(tmp_prefixed_binary_name, 0o550)  # u,g=-w+x


async def async_subprocess_io_multipart(storage, ocaml_diff_role, result_key, write_stdin_function):
    storage_conf = lithops.Storage().storage_config['aws_s3']
    aiob_session = aiobotocore_sess.get_session()
    async with aiob_session.create_client('s3', storage_conf['region_name'],
                                          aws_secret_access_key=storage_conf['secret_access_key'],
                                          aws_access_key_id=storage_conf['access_key_id']) as multipart_client:

        multipart_writer = MultipartUploader(
            multipart_client, storage.bucket, result_key)
        await multipart_writer.setup()

        async def read_stdout(stdout):
            while True:
                buf = await stdout.read(MultipartUploader.PART_SIZE)
                if not buf:
                    break

                await multipart_writer.write(buf)

            await multipart_writer.flush()

        async def read_stderr(stderr):
            while True:
                buf = await stderr.read(2*1024*1024)
                if not buf:
                    break
                print(buf)

        cmd = ' '.join(args.program_and_args)
        print(f'To execute: \'{cmd}\'')
        p = await asyncio.create_subprocess_shell(cmd,
                                                  stdin=asyncio.subprocess.PIPE,
                                                  stdout=asyncio.subprocess.PIPE,
                                                  stderr=asyncio.subprocess.PIPE,
                                                  env=dict(
                                                      os.environ, DIFF=ocaml_diff_role),
                                                  cwd="/tmp")

        await asyncio.gather(
            read_stderr(p.stderr),
            read_stdout(p.stdout),
            write_stdin_function(p.stdin))

        await p.wait()

        if p.returncode != 0:
            raise RuntimeError(
                f"Return code: {p.returncode}")


def run_reducer(results, args, storage):
    download_binary(storage, args)

    async def write_stdin(stdin):
        for cloudobject_path in results:
            partial_result = storage.get_object(
                storage.bucket, cloudobject_path)
            stdin.write(partial_result)
            await stdin.drain()

        stdin.close()
        await stdin.wait_closed()

    session_id = os.environ.get('__LITHOPS_SESSION_ID', '')
    result_key = f'results/{session_id}'

    asyncio.run(async_subprocess_io_multipart(
        storage, 'reducer', result_key, write_stdin), debug=True)

    for cloudobject_path in results:
        storage.delete_object(storage.bucket, cloudobject_path)

    msg = f'lithops storage get {storage.bucket} {result_key} -o RESULTFILENAME'
    print(msg)
    return msg


def print_chunk_info(obj):
    # obj is a CloudObject
    print(obj)
    print(f'part: {obj.part}')
    print(f'data_byte_range: {obj.data_byte_range}')
    print(f'chunk_size: {obj.chunk_size}')


def run_worker(obj, storage, args):
    print_chunk_info(obj)
    download_binary(storage, args)

    async def write_stdin(stdin):
        while True:
            chunk = obj.data_stream.sb.read()
            if not chunk:
                break

            # Discard byte from previous chunk if there is any
            # (this is done so that the program can assume that there is no first byte)
            # if obj.data_stream._plusbytes:
            #    chunk = chunk[1:]

            stdin.write(chunk)
            await stdin.drain()

        stdin.close()
        await stdin.wait_closed()

    session_id = os.environ.get('__LITHOPS_SESSION_ID', '')
    result_key = f'intermediate-results/{session_id}'

    asyncio.run(async_subprocess_io_multipart(
        storage, 'worker', result_key, write_stdin), debug=True)

    return result_key


if __name__ == "__main__":

    ### Parse CLI arguments ###
    parser = argparse.ArgumentParser(
        description='Example: python LaunchJob.py bucketname foldername inputfilepath -- ./MergeOrPair.native merge --lines-per-block 20000')
    parser.add_argument(
        'bucket_name', help='name of the bucket in cloud storage')
    parser.add_argument('input_file_path',
                        help='input file path (local file)')
    parser.add_argument('program_and_args', nargs='*',
                        help='the program to execute without input and output args: ./WhatEver.native --lines-per-block 1000')
    parser.add_argument(
        '--memory', dest='runtime_memory', help='runtime memory in MegaBytes: chunk size is one eigth (1/8) of it', type=int, default=1770)
    args = parser.parse_args()

    # Upload the program first
    if not lithops.Storage().upload_file(args.program_and_args[0], lithops.Storage().bucket,
                                         STORAGE_PROGRAMS_FOLDER + args.program_and_args[0][2:]):
        raise RuntimeError(
            f'Could not upload the program to \'{STORAGE_PROGRAMS_FOLDER}\' subfolder')

    chunk_size_mb = int(args.runtime_memory // 8)
    if chunk_size_mb < 10:
        raise argparse.ArgumentError(
            "Set a chunk size bigger than 10 MB if possible!")

    args.program_and_args.append(f'--chunk-size {chunk_size_mb}')

    ### Launch workers and reducer ###
    fexec_map = lithops.FunctionExecutor(
        log_level='DEBUG', runtime='aiobotocore-vanilla-pyv39:v0')

    inputs_path = f'{args.bucket_name}/{args.input_file_path}'

    map_futures = fexec_map.map(run_worker, [inputs_path],
                                extra_args=[args],
                                runtime_memory=args.runtime_memory,
                                obj_chunk_size=chunk_size_mb*1024*1024)

    intermediate_cloudobjs = fexec_map.get_result(fs=map_futures)
    fexec_map.plot()

    avg_time_sec = sum([x.stats['worker_exec_time']
                       for x in map_futures]) / len(map_futures)
    file_input_size_aprox = chunk_size_mb*len(map_futures)
    print(f'Throughput: { file_input_size_aprox / avg_time_sec : .2f} MB/s')

    # Batch reducer
    fexec_reducer = lithops.FunctionExecutor(
        backend="aws_batch", log_level='DEBUG', runtime='batch_aiobotocore:v0')

    # Lambda reducer
    # fexec_reducer = lithops.FunctionExecutor(
    #    log_level='DEBUG', runtime='aiobotocore-vanilla-pyv39:v0')

    # Increase runtime memory so that we get more vCpus
    reducer_memory = args.runtime_memory * 2

    fexec_reducer.call_async(
        run_reducer, {"results": intermediate_cloudobjs, "args": args}, runtime_memory=reducer_memory)

    print(fexec_reducer.get_result())

    fexec_map.clean()
    fexec_reducer.clean()
