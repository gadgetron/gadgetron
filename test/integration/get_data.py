import os
import sys
import urllib
import hashlib

def md5sum(filename, blocksize=64*1024):
    hsh = hashlib.md5()
    with open(filename, "r+b") as f:
        buf = f.read(blocksize)
        while len(buf) > 0:
            hsh.update(buf)
            buf = f.read(blocksize)
    return hsh.hexdigest()

host = 'http://gadgetrontestdata.s3-website-us-east-1.amazonaws.com'
data = {
    'gre_3d': {
        'simple_gre_out_3d.h5': '60004b2677b9780af9c633e47e7acd1a',
        'meas_MID248_gre_FID30644.dat': '39ac16864627691cf7d84aa2ce13c1ae'
        },
    'radial_phantom': {
        'fixed_radial_mode1.h5':'3a8e10388a3a11683c7e611537c1bd44',
        'golden_radial_mode2.h5':'42c60da3121fa50b8c04e1711b8f4659',
        'meas_MID00133_FID20080_CV_Radial_Fixed_Angle_128_x8_32phs.dat':'58f8de6b6e755c4d3dcd0c7dace3b8f6',
        'meas_MID00135_FID20082_CV_Radial_Golden_Angle_128_512_views.dat':'0326afbb168982f4144704781a08b3ec'
        },
    'rtgrappa': {
        'acc_data_with_device_2.dat':'ac0b59c6c8989c94738e41e2c4b5ec13',
        'grappa_rate2_out.h5':'3d8fa58d851277f8a41b20b02cf82c87'
        },
    'simple_gre': {
        'meas_MiniGadgetron_GRE.dat':'7c5c255522e42367546b4045560afcf8',
        'simple_gre_out.h5':'624ac3178e15e27e52489f330b3fffa5'
        },
    'spiral': {
        'simple_spiral_out.h5':'44be83612c69f008ee71a47fecd3c7ed',
        'meas_MID1132_MiniIRT_spiral_16int_tr500_acc1_10reps_FID13142.dat':'763baf3d7d0acff185ec9a3c85d5a3f3'
        }
}

def main():
    testdir = os.path.dirname(os.path.realpath(sys.argv[0]))
    datadir = os.path.join(testdir, 'data')
    if not os.path.isdir(datadir):
        os.mkdir(datadir)

    print("Storing test data in %s" % datadir)

    for dataset,files in data.items():
        dset_dir = os.path.join(datadir, dataset)
        if not os.path.isdir(dset_dir):
            os.mkdir(dset_dir)
        for fname,checksum in files.items():
            url = '%s/%s/%s' % (host, dataset, fname)
            path = os.path.join(dset_dir, fname)
            tag = '%s/%s' % (dataset, fname)
            localsum = ''
            if os.path.isfile(path):
                print("Computing md5sum of %s..." % tag)
                localsum = md5sum(path)
            if not os.path.isfile(path) or localsum != checksum:
                print("Downloading: %s..." % tag)
                urllib.urlretrieve(url, path)
                print("Verifying: %s..." % tag)
                if md5sum(path) != checksum:
                    print("FAILED!: %s" % tag)

if __name__ == '__main__':
    main()
