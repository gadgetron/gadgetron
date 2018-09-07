#!/usr/bin/python3

import os
import os.path

import sys
import time
import json
import hashlib
import argparse

import urllib.request


class Progress:

    def __init__(self):
        self.start = time.time()
        self.end = None

    def notify(self, blocks, block_size, total_size):

        current_size = blocks * block_size

        if total_size <= current_size:
            self.end = time.time()
            duration = self.end - self.start
            print(' ' * 96, end='\r')
            print("Downloaded {:n} bytes at {:n} bytes per second.".format(current_size, int(total_size / duration)))
        else:
            print(' ' * 96, end='\r')
            print("\t{:n} of {:n} bytes [{:.2%}]".format(current_size, total_size, current_size / total_size), end='\r')


def is_valid(file, digest):

    if not os.path.isfile(file):
        return False

    md5 = hashlib.new('md5')

    with open(file, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            md5.update(chunk)

    return digest == md5.hexdigest()


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

    for entry in entries:

        url = "{}{}".format(args.host, entry['file'])
        destination = os.path.join(args.destination, entry['file'])

        if is_valid(destination, entry['md5']):
            print("Verified: {}".format(destination))
            continue

        print("Downloading file: {}".format(destination))

        os.makedirs(os.path.dirname(destination), exist_ok=True)
        urllib.request.urlretrieve(url, destination, reporthook=Progress().notify)

        if not is_valid(destination, entry['md5']):
            print("Downloaded file {} failed validation.".format(destination))
            sys.exit(1)

    sys.exit(0)


if __name__ == '__main__':
    main()
