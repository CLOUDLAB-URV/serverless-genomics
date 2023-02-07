import lithops
import boto3
import sys

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Please provide a bucket name.")
        print("Usage: python3 abort_mpu.py <bucket_name>")
        sys.exit(1)
        
    bucket_name = sys.argv[1]
    print("Bucket: ", bucket_name)
    
    s3 = lithops.Storage().get_client()
    
    results = s3.list_multipart_uploads(Bucket=bucket_name)
    mpus = results.get('Uploads')
    
    if mpus is not None:
        print("UNFINISHED MULTIPART UPLOADS:")
        for mpu in mpus:
            print("Deleting: ", mpu['Key'])
            s3.abort_multipart_upload(Bucket=bucket_name, Key=mpu['Key'], UploadId=mpu['UploadId'])
        
    print("FINISHED!")