import boto3
import os
from datetime import datetime
from botocore.exceptions import ClientError
from tqdm import tqdm
import threading

class ProgressPercentage(object):
    def __init__(self, filename):
        self._filename = filename
        self._size = float(os.path.getsize(filename))
        self._seen_so_far = 0
        self._lock = threading.Lock()
        self._pbar = tqdm(total=self._size, unit='B', unit_scale=True, desc=filename)

    def __call__(self, bytes_amount):
        # To simplify we'll assume this is hooked up
        # to a single filename for now
        with self._lock:
            self._seen_so_far += bytes_amount
            self._pbar.update(bytes_amount)

    def close(self):
        self._pbar.close()


def check_file_exists_s3(object_name):
    s3_client = boto3.client('s3')
    bucket_name = os.getenv('AWS_BUCKET_NAME')
    try:
        # Attempt to retrieve the metadata of the object
        s3_client.head_object(Bucket=bucket_name, Key=object_name)
        print(f"File '{object_name}' exists in bucket '{bucket_name}'.")
        return True
    except ClientError as e:
        # Check if the error was a 404 (Not Found)
        if e.response['Error']['Code'] == '404':
            print(f"File '{object_name}' does not exist in bucket '{bucket_name}'.")
            return False
        else:
            # Something else went wrong (e.g., permissions issue)
            raise


def upload_file_to_s3(file_name, destination_name=None):
    # If S3 object_name is not specified, use file_name
    bucket_name = os.getenv('AWS_BUCKET_NAME')
    if destination_name is None:
        destination_name = os.path.basename(file_name)

    s3_uri = f"s3://{bucket_name}/{destination_name}"

    if check_file_exists_s3(destination_name):
        return s3_uri

    # Create an S3 client
    s3_client = boto3.client('s3')

    try:
        # Upload the file
        progress = ProgressPercentage(file_name)
        s3_client.upload_file(Filename=file_name,
                              Bucket=bucket_name,
                              Key=destination_name,
                              Callback=progress)
        print(f'File {file_name} uploaded to {destination_name}')

        return s3_uri
    except Exception as e:
        raise Exception(f'An error occurred trying to upload file {file_name} to the AWS storage', e)
    finally:
        progress.close()



def get_today_date():
    today = datetime.today()
    return today.strftime('%Y-%m-%d')
