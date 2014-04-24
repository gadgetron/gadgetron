import os
import sys
import urllib
import hashlib

DATAFILE = "data.txt"
DATADIR = "data"
HOST = 'http://gadgetrontestdata.s3-website-us-east-1.amazonaws.com'

def md5sum(filename, blocksize=64*1024):
    hsh = hashlib.md5()
    with open(filename, "r+b") as f:
        buf = f.read(blocksize)
        while len(buf) > 0:
            hsh.update(buf)
            buf = f.read(blocksize)
    return hsh.hexdigest()

def load_data_file(datafile):
    checksums = {}
    with open(datafile) as f:
        for line in f:
            filepath, checksum = line.split(':')
            checksums[filepath.strip()] = checksum.strip()
    return checksums

def download(url, dest):
    fname = os.path.basename(dest)
    print("Downloading: %s..." % fname)
    urllib.urlretrieve(url, dest)

def main():
    testdir = os.path.dirname(os.path.realpath(sys.argv[0]))
    datadir = os.path.join(testdir, DATADIR)
    datafile = os.path.join(testdir, DATAFILE)
    if not os.path.isdir(datadir):
        os.mkdir(datadir)

    print("Reading list of data from %s" % datafile)
    checksums = load_data_file(datafile)

    print("Storing test data in %s" % datadir)

    # if file exists, verify it. if that fails, download it
    # if file does not exist, download it.
    for dataname,checksum in checksums.items():
        datapath = os.path.join(datadir, dataname)
        parent = os.path.dirname(datapath)
        url = '%s/%s' % (HOST, dataname)
        if os.path.isfile(datapath):
            print("Verifying: %s..." % dataname)
            if md5sum(datapath) != checksum:
                download(url, datapath)
        else:
            if not os.path.isdir(parent):
                os.makedirs(parent)
            download(url, datapath)

if __name__ == '__main__':
    main()
