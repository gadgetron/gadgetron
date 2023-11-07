#!/usr/bin/python3

import pytest
import os
import socket
import hashlib
import subprocess
import yaml
import time

import h5py
import ismrmrd
import numpy


import urllib.error
import urllib.request

from pathlib import Path


# Importing h5py on windows will mess with your environment. When we pass the messed up environment to gadgetron
# child processes, they won't load properly. We're saving our environment here to spare our children from the
# crimes of h5py.
environment = dict(os.environ)


# http://gadgetrondata.blob.core.windows.net/gadgetrontestdata/

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
def fetch_data_file(cache_path: Path, data_host_url:str, cache_disable:bool):
    created_files = []
    
    def _fetch_data_file(filename:str, fileHash:str):
        nonlocal created_files
        url = "{}{}".format(data_host_url, filename)
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
    
    yield _fetch_data_file

    if cache_disable:
        for file in created_files:
            os.remove(file)
            if len(os.listdir(os.path.dirname(file))) == 0:
                os.removedirs(os.path.dirname(file))

@pytest.fixture
def siemens_to_ismrmrd(tmp_path: Path, fetch_data_file):
    def _siemens_to_ismrmrd(fileConfig:dict, section:str):
        input = fetch_data_file(fileConfig['data_file'], fileConfig['data_file_hash'])
        output = os.path.join(tmp_path, os.path.basename(input) + "_" + section + ".h5")

        command = ["siemens_to_ismrmrd", "-X",
                "-f", input, 
                "-m", fileConfig.get('parameter_xml', 'IsmrmrdParameterMap_Siemens.xml'),
                "-x", fileConfig.get('parameter_xsl', 'IsmrmrdParameterMap_Siemens.xsl'),
                "-o", output,
                "-z", str(fileConfig['measurement']),
                fileConfig.get('data_conversion_flag', '')]

        subprocess.run(command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, 
                    cwd=tmp_path)
        
        return output
    
    return _siemens_to_ismrmrd


@pytest.fixture
def send_to_gadgetron(tmp_path: Path, host_url, port):
    def _send_to_gadgetron(fileConfig:dict, input_file:str, section:str):
        output_file = os.path.join(tmp_path, section + ".output.mrd")
        additional_arguments = fileConfig.get('additional_arguments', '')
        configuration = fileConfig['configuration']
        
        print("Passing data to Gadgetron: {} -> {}".format(input_file, output_file)) 
        command = ["gadgetron_ismrmrd_client",
                    "-a", host_url,
                    "-p", port,
                    "-f", input_file,
                    "-o", output_file,
                    "-G", configuration,
                    "-c", configuration]

        if additional_arguments:
            command = command + additional_arguments.split()

        subprocess.run(command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, 
                    env=environment)

        return output_file

    return _send_to_gadgetron

def wait_for_storage_server(port, proc, retries=50):
    for i in range(retries):
        try:
            urllib.request.urlopen(f"http://localhost:{port}/healthcheck")
            return
        except (urllib.error.URLError, urllib.error.HTTPError) as e:
            if i == retries - 1 or proc.poll() is not None:
                raise RuntimeError("Unable to get a successful response from storage server.") from e
            time.sleep(0.2)

def start_storage_server(*, log, port, storage_folder):
    storage_server_environment = environment.copy()
    storage_server_environment["MRD_STORAGE_SERVER_PORT"] = port
    storage_server_environment["MRD_STORAGE_SERVER_STORAGE_CONNECTION_STRING"] = storage_folder
    storage_server_environment["MRD_STORAGE_SERVER_DATABASE_CONNECTION_STRING"] = storage_folder + "/metadata.db"

    retries = 5
    for i in range(retries):
        print("Starting MRD Storage Server on port", port)
        proc = subprocess.Popen(["mrd-storage-server", "--require-parent-pid", str(os.getpid())],
                                stdout=log,
                                stderr=log,
                                env=storage_server_environment)

        try:
            wait_for_storage_server(port, proc)
            return proc
        except:
            # If the process has exited, it might be because the
            # port was in use. This can be because the previous storage server
            # instance was just killed. So we try again.
            if proc.poll() is not None and i < retries:
                time.sleep(1)
            else:
                proc.kill()
                raise

def start_gadgetron_instance(*, log_stdout, log_stderr, port, storage_address, env=environment):
    print("Starting Gadgetron instance on port", port)
    proc = subprocess.Popen(["gadgetron", "-p", port, "-E", storage_address],
                            stdout=log_stdout,
                            stderr=log_stderr,
                            env=env)
    return proc


@pytest.fixture(scope="module", autouse="true")
def start_gadgetron(host_url, port, external, storage_port, tmp_path_factory):
    storage_address = "http://localhost:" + storage_port
    log_path = tmp_path_factory.mktemp("logs")
    storage_path = tmp_path_factory.mktemp("storage")

    if not external:
        print("Starting storage server on port", storage_port)
        with open(os.path.join(log_path, 'storage.log'), 'w') as log:
            storage = start_storage_server(log=log,
                        port=str(storage_port),
                        storage_folder=str(storage_path))

        print("Starting Gadgetron instance on port {} with logs in {}".format(port, log_path))

        with open(os.path.join(log_path, 'gadgetron.log.out'), 'w') as log_stdout:
            with open(os.path.join(log_path, 'gadgetron.log.err'), 'w') as log_stderr:
                instance = start_gadgetron_instance(log_stdout=log_stdout, log_stderr=log_stderr, port=port, storage_address=storage_address)

                yield

                instance.kill()
                storage.kill()

    else:
        yield

Failure = "Failure", 1

def validate_output(*, output_file, reference_file, output_group, reference_group, value_threshold, scale_threshold):
    try:
        # The errors produced by h5py are not entirely excellent. We spend some code here to clear them up a bit.
        def get_group_data(file, group):
            with h5py.File(file, mode='r') as f:
                try:
                    group = group + '/data'
                    return numpy.squeeze(f[group])
                except KeyError:
                    raise RuntimeError("Did not find group '{}' in file {}".format(group, file))

        output_data = get_group_data(output_file, output_group)
        reference_data = get_group_data(reference_file, reference_group)
    except OSError as e:
        return Failure, str(e)
    except RuntimeError as e:
        return Failure, str(e)

    output = output_data[...].flatten().astype('float32')
    reference = reference_data[...].flatten().astype('float32')

    norm_diff = numpy.linalg.norm(output - reference) / numpy.linalg.norm(reference)
    scale = numpy.dot(output, output) / numpy.dot(output, reference)

    if value_threshold < norm_diff:
        return Failure, "Comparing values, norm diff: {} (threshold: {})".format(norm_diff, value_threshold)

    if value_threshold < abs(1 - scale):
        return Failure, "Comparing image scales, ratio: {} ({}) (threshold: {})".format(scale, abs(1 - scale),
                                                                                        scale_threshold)

    return None, "Norm: {:.1e} [{}] Scale: {:.1e} [{}]".format(norm_diff, value_threshold, abs(1 - scale),
                                                                scale_threshold)



# @pytest.mark.parametrize('val', "test_val", ids="test_id")
def test_123(fetch_data_file, siemens_to_ismrmrd, send_to_gadgetron):

    with open("test_cases.yml", 'r') as file:
        cases = yaml.safe_load(file)

    config = cases['cases'][0]

    reconstruction_file = siemens_to_ismrmrd(config['reconstruction'], "reconstruction")
    output_file = send_to_gadgetron(config['reconstruction'], reconstruction_file, "reconstruction")

    for output_image in config['validation']['images']:        
        output_image_data = config['validation']['images'][output_image]
        reference_file = fetch_data_file(output_image_data['reference_file'], output_image_data['reference_file_hash'])
        result, reason = validate_output(output_file=output_file, 
                                        reference_file=reference_file, 
                                        output_group=output_image,
                                        reference_group=output_image_data['reference_image'], 
                                        value_threshold=output_image_data['value_comparison_threshold'], 
                                        scale_threshold=output_image_data['scale_comparison_threshold'])
        
        assert result != Failure, reason

