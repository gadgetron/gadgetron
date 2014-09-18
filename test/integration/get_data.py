import os
import sys
import urllib2
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

def load_checksums(datafile):
    checksums = {}
    with open(datafile) as f:
        for line in f:
            filepath, checksum = line.split(':')
            checksums[filepath.strip()] = checksum.strip()
    return checksums

def download(url, dest):
    furl = urllib2.urlopen(url)
    with open(dest, 'wb') as fdest:
        fdest.write(furl.read())

def main():
    # determine test dir from full path to this script
    testdir = os.path.dirname(os.path.realpath(sys.argv[0]))
    datadir = os.path.join(testdir, DATADIR)
    datafile = os.path.join(testdir, DATAFILE)
    if not os.path.isdir(datadir):
        os.mkdir(datadir)

    print("Reading list of data from %s" % datafile)
    try:
        checksums = load_checksums(datafile)
    except IOError:
        print("Failed to read %s" % datafile)
        return

    print("Storing test data in %s" % datadir)

    for dataname,checksum in checksums.items():
        datapath = os.path.join(datadir, dataname)
        parent = os.path.dirname(datapath)
        if not os.path.isdir(parent):
            os.makedirs(parent)
        url = '%s/%s' % (HOST, dataname)

        print("Verifying: %s..." % dataname)
        # if file is missing or its checksum doesn't match, download it
        if not os.path.isfile(datapath) or md5sum(datapath) != checksum:
            print("Downloading: %s..." % dataname)
            try:
                download(url, datapath)
            except urllib2.HTTPError, e:
                print("HTTP Error: %d %s" % (e.code, url))
            except urllib2.URLError, e:
                print("URL Error: %s - %s" % (e.reason, url))

if __name__ == '__main__':
    main()
