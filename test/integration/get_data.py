#!/usr/bin/python3

import os
import os.path

import sys
import time
import json
import hashlib
import argparse
import functools

import urllib.error
import urllib.request

from concurrent.futures import ThreadPoolExecutor

def is_valid(file, digest):

    if not os.path.isfile(file):
        return False

    md5 = hashlib.new('md5')

    with open(file, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            md5.update(chunk)

    return digest == md5.hexdigest()

def urlretrieve(url, filename, retries=3):
    if retries <= 0:
        raise RuntimeError("Download from {} failed".format(url))
    try:
        with urllib.request.urlopen(url, timeout=5) as connection:
            with open(filename,'wb') as f:
                for chunk in iter(lambda : connection.read(1024*1024), b''):
                    f.write(chunk)
    except urllib.error.URLError:
        urlretrieve(url, filename, retries=retries-1)





def main():

    parser = argparse.ArgumentParser(description="Gadgetron Integration Test Data Download Script",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-l', '--list', type=str, default='data.json',
                        help="List of data files to download.")

    parser.add_argument('-d', '--destination', type=str, default='data',
                        help="Folder in which to write downloaded data.")

    parser.add_argument('-H', '--host', default='http://gadgetrondata.blob.core.windows.net/gadgetrontestdata/',
                        help="Host from which to download the data.")

    args = parser.parse_args()

    with open(args.list, 'r') as list:
        entries = json.load(list)

    def download_entry(entry):
        url = "{}{}".format(args.host, entry['file'])
        destination = os.path.join(args.destination, entry['file'])

        if is_valid(destination, entry['md5']):
            print("Verified: {}".format(destination))
            return 

        print("Downloading file: {}".format(url))

        os.makedirs(os.path.dirname(destination), exist_ok=True)
        urlretrieve(url,destination)

        if not is_valid(destination, entry['md5']):
            raise(RuntimeError("Downloaded file {} failed validation.".format(destination)))

        print("File saved as: {}".format(destination))
            
    with ThreadPoolExecutor() as executor:
        executor.map(download_entry,entries)    

if __name__ == '__main__':
    main()
