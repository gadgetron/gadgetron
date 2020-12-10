#!/usr/bin/python3

import os
import os.path

import sys
import time
import json
import hashlib
import argparse
import functools

import urllib.request
import requests
from multiprocessing import Pool

from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

retry_strategy = Retry(total=3)
adapter = HTTPAdapter(max_retries=retry_strategy)
http = requests.Session()
http.mount("https://", adapter)
http.mount("http://", adapter)

def is_valid(file, digest):
    if not os.path.isfile(file):
        return False

    md5 = hashlib.new('md5')

    with open(file, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            md5.update(chunk)

    return digest == md5.hexdigest()

def download_file(entry):
    url = entry['url']
    destination = entry['destination']
    os.makedirs(os.path.dirname(destination), exist_ok=True)

    if is_valid(destination, entry['md5']):
        print("Verified: {}".format(destination))
        return

    print("Downloading file: {}".format(destination))
 
    r = http.get(url, stream=True, timeout=10)
    if r.status_code == 200:
        with open(destination, 'wb') as f:
            for chunk in r:
                f.write(chunk)

    if not is_valid(destination, entry['md5']):
        print("Downloaded file {} failed validation.".format(destination))
        sys.exit(1)


def __update_handler(current_size, total_size):
    print(' ' * 64, end='\r')
    print("\t{:n} of {:n} bytes [{:.2%}]".format(current_size, total_size, current_size / total_size), end='\r')


def __finish_handler(total_size, duration):
    print(' ' * 64, end='\r')
    print("Downloaded {:n} bytes at {:n} bytes per second.".format(total_size, int(total_size / duration)))


def __progress_notification_handler(on_update, on_finish, start, blocks, block_size, total_size):
    current_size = blocks * block_size

    if total_size <= current_size:
        on_finish(total_size, time.time() - start)
    else:
        on_update(current_size, total_size)



def main():

    parser = argparse.ArgumentParser(description="Gadgetron Integration Test Data Download Script",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-l', '--list', type=str, default='data.json',
                        help="List of data files to download.")

    parser.add_argument('-d', '--destination', type=str, default='data',
                        help="Folder in which to write downloaded data.")

    parser.add_argument('-H', '--host', default='http://gadgetrondata.blob.core.windows.net/gadgetrontestdata/',
                        help="Host from which to download the data.")

    parser.add_argument('--mute-download-progress', dest='progress', action='store_const',
                        const=functools.partial(__progress_notification_handler, lambda *_: None, __finish_handler),
                        default=functools.partial(__progress_notification_handler, __update_handler, __finish_handler),
                        help="Mute download progress messages.")

    args = parser.parse_args()

    with open(args.list, 'r') as list:
        entries = json.load(list)

    resolved_entries = []
    for entry in entries:
        url = "{}{}".format(args.host, entry['file'])
        destination = os.path.join(args.destination, entry['file'])
        resolved_entries.append({ 'url': url, 'destination': destination, 'md5': entry['md5'] })

    print('Starting threadpool')
    with Pool(30) as p:
        p.map(download_file, resolved_entries)

    #ThreadPool(3).imap_unordered(download_file, entries)

    sys.exit(0)


if __name__ == '__main__':
    main()
