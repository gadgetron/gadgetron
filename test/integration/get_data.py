#!/usr/bin/python3

import os
import os.path

import sys
import time
import json
import hashlib
import functools
import argparse

import requests
from requests.adapters import TimeoutSauce
import backoff
from multiprocessing import Pool

REQUESTS_TIMEOUT_SECONDS = 10
REQUESTS_MAX_RETRIES = 3

class CustomTimeout(TimeoutSauce):
    def __init__(self, *args, **kwargs):
        if kwargs["connect"] is None:
            kwargs["connect"] = REQUESTS_TIMEOUT_SECONDS
        if kwargs["read"] is None:
            kwargs["read"] = REQUESTS_TIMEOUT_SECONDS
        super().__init__(*args, **kwargs)

# Setting the global time out
requests.adapters.TimeoutSauce = CustomTimeout

# Retry decorator with exponential wait.
retry_timeout = backoff.on_exception(
    wait_gen=backoff.expo,
    exception=(
        requests.exceptions.Timeout,
        requests.exceptions.ConnectionError
    ),
    max_tries=REQUESTS_MAX_RETRIES,
)

def is_valid(file, digest):
    if not os.path.isfile(file):
        return False

    md5 = hashlib.new('md5')

    with open(file, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            md5.update(chunk)

    return digest == md5.hexdigest()

@retry_timeout
def download_file(entry):
    url = entry['url']
    destination = entry['destination']
    os.makedirs(os.path.dirname(destination), exist_ok=True)

    if is_valid(destination, entry['md5']):
        print("Verified: {}".format(destination))
        return

    print("Downloading file: {}".format(destination))
 
    r = requests.get(url, stream=True)
    if r.status_code == 200:
        with open(destination, 'wb') as f:
            for chunk in r:
                f.write(chunk)

    if not is_valid(destination, entry['md5']):
        print("Downloaded file {} failed validation.".format(destination))
        sys.exit(1)


def main():

    parser = argparse.ArgumentParser(description="Gadgetron Integration Test Data Download Script",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-l', '--list', type=str, default='data.json',
                        help="List of data files to download.")

    parser.add_argument('-d', '--destination', type=str, default='data',
                        help="Folder in which to write downloaded data.")

    parser.add_argument('-H', '--host', default='http://gadgetrondata.blob.core.windows.net/gadgetrontestdata/',
                        help="Host from which to download the data.")

    parser.add_argument('-t', '--threads', type=int, default=50,
                        help="Number of download threads")

    args = parser.parse_args()

    with open(args.list, 'r') as list:
        entries = json.load(list)

    resolved_entries = []
    for entry in entries:
        url = "{}{}".format(args.host, entry['file'])
        destination = os.path.join(args.destination, entry['file'])
        resolved_entries.append({ 'url': url, 'destination': destination, 'md5': entry['md5'] })

    with Pool(args.threads) as p:
        p.map(download_file, resolved_entries)

    sys.exit(0)


if __name__ == '__main__':
    main()
