#!/usr/bin/python3

import pytest
import os
import socket
import hashlib
import subprocess

import urllib.error
import urllib.request

from pathlib import Path

# http://gadgetrondata.blob.core.windows.net/gadgetrontestdata/

@pytest.fixture
def data_dir():
    current_dir = Path(os.path.dirname(__file__))
    return current_dir/"data"

def loaddata():
    return (["cmr/CineBinning/meas_MID838_PK_rt_test_2slice_FID22519/meas_MID838_PK_rt_test_2slice_FID22519.dat"], ["make1"])

def loaddataX():
    return [("test1", "make1"), ("test2", "make3"), ("test3", "make5")]


def loaddata2():
    global test_val, test_id, test_g
    test_val, test_id = loaddata()
    test_g = loaddataX()

loaddata2()

@pytest.fixture(params=test_val)
def named(request):
    return request.param

def calc_mdf5(file):
    md5 = hashlib.new('md5')
    
    with open(file, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            md5.update(chunk)
    return md5.hexdigest()

def is_valid(file, digest):

    if not os.path.isfile(file):
        return False

    return digest == calc_mdf5(file) 

def urlretrieve(url, filename, retries=5):
    if retries <= 0:
        raise RuntimeError("Download from {} failed".format(url))
    try:
        with urllib.request.urlopen(url, timeout=60) as connection:                        
            with open(filename,'wb') as f:
                for chunk in iter(lambda : connection.read(1024*1024), b''):
                    f.write(chunk)
            return connection.headers["Content-MD5"]
    except (urllib.error.URLError, ConnectionResetError, socket.timeout) as exc:
        print("Retrying connection for file {}, reason: {}".format(filename, str(exc)))
        urlretrieve(url, filename, retries=retries-1)

@pytest.fixture
def convertFile(tmp_path: Path):
    def _convertFile(input:str, fileConfig:dict):
        working_dir = os.path.join(tmp_path, os.path.basename(input) + "output")
        os.makedirs(working_dir, exist_ok=True)
        output = os.path.join(working_dir, os.path.basename(input) + ".h5")

        command = ["siemens_to_ismrmrd", "-X",
                "-f", input, 
                "-m", fileConfig.get('parameter_xml', 'IsmrmrdParameterMap_Siemens.xml'),
                "-x", fileConfig.get('parameter_xsl', 'IsmrmrdParameterMap_Siemens.xsl'),
                "-o", output,
                "-z", fileConfig['measurement'],
                fileConfig.get('data_conversion_flag', '')]

        print(command)
        subprocess.run(command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, 
                    cwd=working_dir)
        
        return output
    
    return _convertFile

@pytest.fixture
def dataFile(cache_path: Path, host_url:str, cache_disable:bool):
    created_files = []
    
    def _dataFile(filename:str, fileHash:str):
        nonlocal created_files
        url = "{}{}".format(host_url, filename)
        destination = os.path.join(cache_path, filename)

        if not os.path.exists(destination):
            print("Downloading file: {}".format(url)) 

            os.makedirs(os.path.dirname(destination), exist_ok=True)
            urlretrieve(url, destination)
            
        else:
            print("File already exists: {}".format(url)) 

        if not is_valid(destination, fileHash):
            raise(RuntimeError("Downloaded file {} failed validation. Expected MD5 {}. Actual MD5 {}".format(destination, fileHash, calc_mdf5(destination))))

        created_files.append(destination)

        return destination
    
    yield _dataFile

    if cache_disable:
        for file in created_files:
            os.remove(file)
            if len(os.listdir(os.path.dirname(file))) == 0:
                os.removedirs(os.path.dirname(file))






# def tearDownModule():
#     raise Exception(f"Run failed")



def localwork(data_dir: str):
    print(data_dir)

@pytest.mark.parametrize('val', test_val, ids=test_id)
def test_123(val: str, dataFile, convertFile):
    path1 = dataFile(val, "7f52ef8abeb44b7a23648b195307ccc8")
    path2 = dataFile("cmr/CineBinning/meas_MID838_PK_rt_test_2slice_FID22519/cmr_cine_binning_2slice_ref_20220817.mrd", "ddd08991b7837dc37d7457cc74fb40cb")
    
    config = {'measurement':'2', 'data_conversion_flag':'--flashPatRef'}

    path3 = convertFile(path1, config)
    print(path1, path2, path3)
    assert 1
