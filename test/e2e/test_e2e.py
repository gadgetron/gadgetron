#!/usr/bin/python3

import pytest
import os
import socket
import hashlib
import subprocess
import yaml
import time
import string
import re
import sys
import itertools
import io

import h5py
import ismrmrd
import numpy
import glob
import json

import urllib.error
import urllib.request

from pathlib import Path
from typing import Dict, List, Callable


def calc_mdf5(file: Path) -> str:
    md5 = hashlib.new('md5')
    
    with open(file, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            md5.update(chunk)
    return md5.hexdigest()

def is_valid(file: Path, digest: str) -> bool:
    if not os.path.isfile(file):
        return False

    return digest == calc_mdf5(file) 

def urlretrieve(url: str, filename: str, retries: int = 5) -> str:
    if retries <= 0:
        pytest.fail(RuntimeError("Download from {} failed".format(url)))
    try:
        with urllib.request.urlopen(url, timeout=60) as connection:                        
            with open(filename,'wb') as f:
                for chunk in iter(lambda : connection.read(1024*1024), b''):
                    f.write(chunk)
            return connection.headers["Content-MD5"]
    except (urllib.error.URLError, ConnectionResetError, socket.timeout) as exc:
        print("Retrying connection for file {}, reason: {}".format(filename, str(exc)))
        return urlretrieve(url, filename, retries=retries-1)

@pytest.fixture
def fetch_data_file(cache_path: Path, data_host_url: str, cache_disable: bool, tmp_path: Path) -> Callable:
    created_files = []
    
    def _fetch_data_file(filename: str, fileHash: str) -> str:
        nonlocal created_files
        url = "{}{}".format(data_host_url, filename)

        if cache_disable:
            destination = os.path.join(tmp_path, filename)
        else:
            destination = os.path.join(cache_path, filename)

        if not os.path.exists(destination):
            print("Downloading file: {}".format(url)) 

            os.makedirs(os.path.dirname(destination), exist_ok=True)
            urlretrieve(url, destination)
            
        else:
            print("File already exists: {}".format(destination)) 

            if not is_valid(destination, fileHash):
                print("Downloaded file {} MD5 mismatch, forcing download. Expected MD5 {}. Actual MD5 {}".format(destination, fileHash, calc_mdf5(destination)))

            os.makedirs(os.path.dirname(destination), exist_ok=True)
            urlretrieve(url, destination)

        if not is_valid(destination, fileHash):
            pytest.fail("Downloaded file {} failed validation. Expected MD5 {}. Actual MD5 {}".format(destination, fileHash, calc_mdf5(destination)))

        created_files.append(destination)

        return destination
    
    yield _fetch_data_file

    if cache_disable:
        for file in created_files:
            if os.path.isfile(file):
                os.remove(file)
                if len(os.listdir(os.path.dirname(file))) == 0:
                    os.removedirs(os.path.dirname(file))

@pytest.fixture
def siemens_to_ismrmrd(tmp_path: Path, fetch_data_file: Callable) -> Callable:
    def _siemens_to_ismrmrd(fileConfig: Dict[str, str], section: str) -> str:
        if 'copy_file' in fileConfig:
            output = fetch_data_file(fileConfig['copy_file'], fileConfig['copy_file_hash'])
            print("Using output file {}".format(output)) 

        else:
            input = fetch_data_file(fileConfig['data_file'], fileConfig['data_file_hash'])
            output = os.path.join(tmp_path, os.path.basename(input) + "_" + section + ".h5")

            print("Convert file {} to {}".format(input, output)) 

            command = ["siemens_to_ismrmrd", "-X",
                    "-f", input, 
                    "-m", fileConfig.get('parameter_xml', 'IsmrmrdParameterMap_Siemens.xml'),
                    "-x", fileConfig.get('parameter_xsl', 'IsmrmrdParameterMap_Siemens.xsl'),
                    "-o", output,
                    "-z", str(fileConfig['measurement']),
                    fileConfig.get('data_conversion_flag', '')]

            with open(os.path.join(tmp_path, os.path.basename(input) + "_" + section + '.log.out'), 'w') as log_stdout:
                with open(os.path.join(tmp_path, os.path.basename(input) + "_" + section + '.log.err'), 'w') as log_stderr:
                    result = subprocess.run(command,
                                stdout=log_stdout,
                                stderr=log_stderr, 
                                cwd=tmp_path)
            
            if result.returncode != 0:
                pytest.fail("siemens_to_ismrmrd failed with return code {}".format(result.returncode))

        return output
    
    return _siemens_to_ismrmrd

def query_gadgetron_capabilities(info_string: str) -> Dict[str, str]:
    print("Querying Gadgetron capabilities...")

    value_pattern = r"(?:\s*):(?:\s+)(?P<value>.*)?"

    capability_markers = {
        'version': "Version",
        'build': "Git SHA1",
        'memory': "System Memory size",
        'python': "Python Support",
        'julia': "Julia Support",
        'matlab': "Matlab Support",
        'cuda': "CUDA Support"
    }

    plural_capability_markers = {
        'cuda_memory': "Total amount of global GPU memory",
        'cuda_devices': "Number of CUDA capable devices"
    }

    def find_value(marker):
        pattern = re.compile(marker + value_pattern, re.IGNORECASE)
        match = pattern.search(info_string)

        if not match:
            pytest.fail("Failed to parse Gadgetron info string; Gadgetron capabilities could not be determined.")

        return match['value']

    def find_plural_values(marker):
        pattern = re.compile(marker + value_pattern, re.IGNORECASE)
        return [match['value'] for match in pattern.finditer(info_string)]

    capabilities = {key: find_value(marker) for key, marker in capability_markers.items()}
    capabilities.update({key: find_plural_values(marker) for key, marker in plural_capability_markers.items()})

    return capabilities

@pytest.fixture
def check_requirements(host_url: str, port: int, external: bool, ignore_requirements: List[str], run_tag: List[str], tmp_path: Path) -> Callable:
    def _check_requirements(requirements: Dict[str, str], tags: List[str], local: bool = False, retries: int = 50) -> None:
        def rules_from_reqs(section):
            class Rule:
                def __init__(self, capability, validator, message):
                    self.capability = capability
                    self.validator = validator
                    self.message = message

                def is_satisfied(self, capabilities):
                    value = capabilities.get(self.capability)
                    return self.validator(value)
        
            def parse_memory(string):
                pattern = re.compile(r"(?P<value>\d+)(?: MB)?")
                match = pattern.search(string)
                return float(match['value'])

            def is_enabled(value):
                return value in ['YES', 'yes', 'True', 'true', '1']

            def has_more_than(target):
                return lambda value: parse_memory(str(target)) <= parse_memory(value)

            def each(validator):
                return lambda values: all([validator(value) for value in values])

            rules = [
                ('matlab_support', lambda req: Rule('matlab', is_enabled, "MATLAB support required.")),
                ('python_support', lambda req: Rule('python', is_enabled, "Python support required.")),
                ('julia_support', lambda req: Rule('julia', is_enabled, "Julia support required.")),
                ('system_memory', lambda req: Rule('memory', has_more_than(req), "Not enough system memory.")),
                ('gpu_support', lambda req: Rule('cuda', is_enabled, "CUDA support required.")),
                ('gpu_support', lambda req: Rule('cuda_devices', each(has_more_than(req)), "Not enough CUDA devices.")),
                ('gpu_memory', lambda req: Rule('cuda_memory', each(has_more_than(req)), "Not enough graphics memory."))
            ]

            return [(key, rule(section[key])) for key, rule in rules if key in section]

        if run_tag != "" and run_tag not in tags:
            pytest.skip("Test missing required tag.")    

        if 'skip' in tags:
            pytest.skip("Test was marked as skipped")

        if local or external: 
            command = ["gadgetron", "--info"]
        else:
            command = ["gadgetron_ismrmrd_client",
                        "-a", host_url,
                        "-p", str(port),
                        "-q", "-Q", "gadgetron::info"]


        sleep_time = 0.2
        for i in range(retries):
            with open(os.path.join(tmp_path, 'gadgetron_info.log'), 'w') as log_stdout:
                info_process = subprocess.run(command, stdout=log_stdout, stderr=log_stdout, universal_newlines=True)            

                if info_process.returncode:
                    if i == retries - 1:
                        pytest.fail("Querying gadgetron info failed with return code {}".format(info_process.returncode))

                    time.sleep(sleep_time)
                    sleep_time += 0.2
                else:
                    break

        with open(os.path.join(tmp_path, 'gadgetron_info.log'), 'w') as log_stdout:
            info_process = subprocess.run(command, stdout=log_stdout, stderr=log_stdout, universal_newlines=True)            
            if info_process.returncode:
                pytest.fail("Querying gadgetron info failed with return code {}".format(info_process.returncode))
            
        with open(os.path.join(tmp_path, 'gadgetron_info.log')) as log_stdout:
            info_string = log_stdout.read()

        # info_string = subprocess.check_output(command, stderr=subprocess.STDOUT, universal_newlines=True)
        capabilities = query_gadgetron_capabilities(info_string)

        reqs = rules_from_reqs(requirements)

        for rule in [rule for key, rule in reqs]:
            # Ignore rules that are in the ignore_requirements list
            if rule.capability in ignore_requirements:
                continue

            if not rule.is_satisfied(capabilities):
                pytest.skip(rule.message)    

    return _check_requirements        

template_path = './config'


def get_config(fileConfig: Dict[str, str], tmp_path: str, section: str) -> (str, List[str]):
    if 'template' in fileConfig:
        template_file = os.path.join(template_path, fileConfig['template'])
        configuration = os.path.join(tmp_path, section + '.config.xml')

        config_option = ["-G", section, 
                        "-C", configuration]

        with open(template_file, 'r') as input:
            with open(configuration, 'w') as output:
                output.write(
                    string.Template(input.read()).substitute(
                        test_folder=os.path.abspath(tmp_path),
                        # Expand substitution list as needed.
                    )
                )

    else:
        configuration = fileConfig['configuration']

        config_option = ["-G", configuration, 
                        "-c", configuration]
        
    return configuration, config_option


@pytest.fixture
def send_to_gadgetron(tmp_path: Path, host_url: str, port: int) -> Callable:
    def _send_to_gadgetron(fileConfig: Dict[str, str], input_file: str, section: str) -> str:
        output_file = os.path.join(tmp_path, section + ".output.mrd")
        additional_arguments = fileConfig.get('additional_arguments', '')

        _, config_option = get_config(fileConfig, tmp_path, section)
        
        print("Passing data to Gadgetron: {} -> {}".format(input_file, output_file)) 
        command = ["gadgetron_ismrmrd_client",
                    "-a", host_url,
                    "-p", str(port),
                    "-f", input_file,
                    "-o", output_file]
        
        command += config_option

        if additional_arguments:
            command = command + additional_arguments.split()

        with open(os.path.join(tmp_path, section + '.client.log'), 'w') as log_stdout:
            result = subprocess.run(command,
                        stdout=log_stdout,
                        stderr=log_stdout)

        if result.returncode != 0:
            pytest.fail("gadgetron_ismrmrd_client failed with return code {}".format(result.returncode))

        return output_file

    return _send_to_gadgetron

@pytest.fixture
def stream_to_gadgetron(tmp_path: Path, storage_port: int, external: bool) -> Callable:
    def _stream_to_gadgetron(fileConfig: Dict[str, str], input_file: str, section: str) -> str:
        output_file = os.path.join(tmp_path, section + ".output.mrd")
        storage_address = "http://localhost:" + str(storage_port)

        input_adapter = fileConfig.get('input_adapter', 'ismrmrd_hdf5_to_stream')
        output_adapter = fileConfig.get('output_adapter', 'ismrmrd_stream_to_hdf5')
        output_group = fileConfig.get('output_group', '')
        streamConfig = fileConfig.get('stream', [])
        
        if not output_group:
            output_group = fileConfig.get('configuration')

        if not streamConfig:
            configuration, _ = get_config(fileConfig, tmp_path, section)
            streamConfig = [ {"configuration": configuration} ]

        stream_command = f"{input_adapter} -i {input_file} --use-stdout"

        for stream in streamConfig:
            stream_command += f" | gadgetron -E {storage_address} --from_stream -c {stream['configuration']} {stream.get('args', '')}"
            
        stream_command += f" | {output_adapter} --use-stdin -o {output_file} -g {output_group}"
        
        # Some stream arguments use ${test_folder} directly so this will provide support for that.
        stream_command = stream_command.replace('${test_folder}', os.path.abspath(tmp_path))

        split_cmd = ['bash', '-c', stream_command]

        print("Streaming data to Gadgetron: {} -> {}".format(input_file, output_file)) 

        with open(os.path.join(tmp_path, section + '_gadgetron.log.out'), 'w') as log_stdout:
            with open(os.path.join(tmp_path, section + '_gadgetron.log.err'), 'w') as log_stderr:
                result = subprocess.run(split_cmd,
                            stdout=log_stdout,
                            stderr=log_stderr)

        if result.returncode != 0:
            pytest.fail("stream processing failed with return code {}".format(result.returncode))

        return output_file

    return _stream_to_gadgetron

def wait_for_storage_server(port: int, proc: subprocess.Popen, retries: int = 50) -> None:
    sleep_time = 0.2
    for i in range(retries):
        try:
            urllib.request.urlopen(f"http://localhost:{port}/healthcheck")
            return
        except (urllib.error.URLError, urllib.error.HTTPError, ConnectionResetError, ConnectionRefusedError, socket.timeout) as e:
            if i == retries - 1 or proc.poll() is not None:
                pytest.fail("Unable to get a successful response from storage server. {}".format(e))
            time.sleep(sleep_time)
            sleep_time += 0.2

def start_storage_server(*, log: str, port: int, storage_folder: str) -> None:
    storage_server_environment = os.environ.copy()
    storage_server_environment["MRD_STORAGE_SERVER_PORT"] = str(port)
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

def start_gadgetron_instance(*, log_stdout: io.TextIOWrapper, log_stderr: io.TextIOWrapper, port: int, storage_address: str, env = None) -> subprocess.Popen:
    print("Starting Gadgetron instance on port", port)
    proc = subprocess.Popen(["gadgetron", "-p", str(port), "-E", str(storage_address)],
                            stdout=log_stdout,
                            stderr=log_stderr,
                            env=env)
    return proc


@pytest.fixture
def start_storage(external: bool, storage_port: int, tmp_path: Path) -> None:
    def _start_storage():
        storage_path = os.path.join(tmp_path, "storage")
        print(storage_path)
        os.mkdir(storage_path)

        storage = None

        if not external:
            print("Starting storage server on port", storage_port)
            with open(os.path.join(tmp_path, 'storage.log'), 'w') as log:
                storage = start_storage_server(log=log,
                            port=int(storage_port),
                            storage_folder=str(storage_path))
                
        return storage 

    return _start_storage

@pytest.fixture
def start_gadgetron_sever(port: int, external: bool, storage_port: int, tmp_path: Path) -> Callable:
    instance = None
    def _start_gadgetron() -> None:
        nonlocal instance

        storage_address = "http://localhost:" + str(storage_port)

        if not external:
            print("Starting Gadgetron instance on port {} with logs in {}".format(port, tmp_path))

            with open(os.path.join(tmp_path, 'gadgetron.log.out'), 'w') as log_stdout:
                with open(os.path.join(tmp_path, 'gadgetron.log.err'), 'w') as log_stderr:
                    instance = start_gadgetron_instance(log_stdout=log_stdout, log_stderr=log_stderr, port=port, storage_address=storage_address)

    yield _start_gadgetron

    if instance != None:
        instance.kill()


@pytest.fixture
def start_gadgetron_sever_with_additional_nodes(tmp_path: Path, port: int, storage_port: int, external: bool) -> Callable:
    instances = []

    def _start_gadgetron_with_additional_nodes(fileConfig: Dict[str, str]) -> None:
        nonlocal instances

        if external:
            return

        base_port = 9050
        number_of_nodes = 2

        if 'distributed' in fileConfig:
            base_port = int(fileConfig['distributed'].get('node_port_base', 9050))
            number_of_nodes = int(fileConfig['distributed'].get('nodes', 2))

        storage_address = "http://localhost:" + str(storage_port)
        env=os.environ.copy()

        ids = range(number_of_nodes)
        print("Will start additional Gadgetron workers on ports:", *map(lambda idx: base_port + idx, ids))

        worker_list = []

        for id in ids:
            instance_port = str(base_port + id)
            with open(os.path.join(tmp_path, 'gadgetron_worker' + instance_port + '.log.out'), 'w') as log_stdout:
                with open(os.path.join(tmp_path, 'gadgetron_worker' + instance_port + '.log.err'), 'w') as log_stderr:
                    instances.append(start_gadgetron_instance(log_stdout=log_stdout, log_stderr=log_stderr, port=instance_port, storage_address=storage_address))
                    worker_list=worker_list + ['localhost:' + str(instance_port)]

        env["GADGETRON_REMOTE_WORKER_COMMAND"] = "echo " + json.dumps(worker_list)

        print("Setting env to", env["GADGETRON_REMOTE_WORKER_COMMAND"])

        print("Starting main Gadgetron instance on port {} with logs in {}".format(port, tmp_path))

        with open(os.path.join(tmp_path, 'gadgetron.log.out'), 'w') as log_stdout:
            with open(os.path.join(tmp_path, 'gadgetron.log.err'), 'w') as log_stderr:
                instances.append(start_gadgetron_instance(log_stdout=log_stdout, log_stderr=log_stderr, port=port, storage_address=storage_address, env=env))

    yield _start_gadgetron_with_additional_nodes

    for instance in instances:
        instance.kill()


def validate_output(*, output_file: str, reference_file: str, output_group: str, reference_group: str, value_threshold: float, scale_threshold: float) -> None:
    try:
        # The errors produced by h5py are not entirely excellent. We spend some code here to clear them up a bit.
        def get_group_data(file, group):
            with h5py.File(file, mode='r') as f:
                try:
                    group = group + '/data'
                    return numpy.squeeze(f[group])
                except KeyError:
                    pytest.fail("Did not find group '{}' in file {}".format(group, file))

        output_data = get_group_data(output_file, output_group)
        reference_data = get_group_data(reference_file, reference_group)
    except OSError as e:
        pytest.fail(str(e))
    except RuntimeError as e:
        pytest.fail(str(e))

    output = output_data[...].flatten().astype('float32')
    reference = reference_data[...].flatten().astype('float32')

    norm_diff = numpy.linalg.norm(output - reference) / numpy.linalg.norm(reference)
    scale = numpy.dot(output, output) / numpy.dot(output, reference)

    if value_threshold < norm_diff:
        pytest.fail("Comparing values, norm diff: {} (threshold: {})".format(norm_diff, value_threshold))

    if value_threshold < abs(1 - scale):
        pytest.fail("Comparing image scales, ratio: {} ({}) (threshold: {})".format(scale, abs(1 - scale),
                                                                                        scale_threshold))

def validate_dataset(*, dataset_file: str, reference_file: str, dataset_group: str, reference_group: str) -> None:
    try:
        dataset_file = ismrmrd.File(dataset_file, 'r')
    except OSError as e:
        pytest.fail("Failed to read dataset file '{}'".format(dataset_file))

    try:
        reference_file = ismrmrd.File(reference_file, 'r')
    except OSError as e:
        pytest.fail("Failed to read reference file '{}'".format(reference_file))

    header = dataset_file[dataset_group].header
    ref_header = reference_file[reference_group].header
    if not dataset_file[dataset_group].header == reference_file[reference_group].header:
        import deepdiff
        diff = deepdiff.diff.DeepDiff(header, ref_header)
        print(diff.pretty())
        pytest.fail("Dataset header did not match reference header")

    for attribute in ['acquisitions', 'waveforms', 'images']:
        dataset = getattr(dataset_file[dataset_group], attribute) or []
        reference = getattr(reference_file[reference_group], attribute) or []

        if not list(dataset) == list(reference):
            pytest.fail("Dataset {attr} did not match reference {attr}".format(attr=attribute))

@pytest.fixture
def validate_dataset_output(tmp_path: Path, fetch_data_file: Callable) -> Callable:
    def _validate_dataset_output(fileConfig: Dict[str, str]) -> None:
        
        dataset_prefix = os.path.join(tmp_path, fileConfig['dataset_prefix'])
        dataset_files = glob.glob(dataset_prefix + "*")

        print(dataset_files)

        if len(dataset_files) == 0:
            pytest.fail("Found no dataset with prefix: {}".format(dataset_prefix))

        if len(dataset_files) > 1:
            pytest.fail("Too many datasets with prefix: {}".format(dataset_prefix))

        reference_file = fetch_data_file(fileConfig['reference_file'], fileConfig['reference_file_hash'])

        validate_dataset(dataset_file=dataset_files[0],
                        dataset_group=fileConfig.get('dataset_group', 'dataset'),
                        reference_file=reference_file,
                        reference_group=fileConfig.get('reference_group', 'dataset'))

    return _validate_dataset_output

@pytest.fixture
def validate_output_data(fetch_data_file: Callable) -> Callable:
    def _validate_output_data(fileConfig: Dict[str, str], output_file: str, output_image: str) -> None:
        reference_file = fetch_data_file(fileConfig['reference_file'], fileConfig['reference_file_hash'])
        validate_output(output_file=output_file, 
                        reference_file=reference_file, 
                        output_group=output_image,
                        reference_group=fileConfig['reference_image'], 
                        value_threshold=fileConfig['value_comparison_threshold'], 
                        scale_threshold=fileConfig['scale_comparison_threshold'])

        if not fileConfig.get('disable_image_header_test', False):
            validate_image_header(output_file=output_file,                         
                                    reference_file=reference_file, 
                                    output_group=output_image,
                                    reference_group=fileConfig['reference_image'])
    
    return _validate_output_data


def validate_image_header(*, output_file: str, reference_file: str, output_group: str, reference_group: str) -> None:
    def equals():
        return lambda out, ref: out == ref

    def approx(threshold=1e-6):
        return lambda out, ref: abs(out - ref) <= threshold

    def ignore():
        return lambda out, ref: True

    def each(rule):
        return lambda out, ref: all(rule(out, ref) for out, ref in itertools.zip_longest(out, ref))

    header_rules = {
        'version': equals(),
        'data_type': equals(),
        'flags': equals(),
        'measurement_uid': equals(),
        'matrix_size': each(equals()),
        'field_of_view': each(approx()),
        'channels': equals(),
        'position': each(approx()),
        'read_dir': each(approx()),
        'phase_dir': each(approx()),
        'slice_dir': each(approx()),
        'patient_table_position': each(approx()),
        'average': equals(),
        'slice': equals(),
        'contrast': equals(),
        'phase': equals(),
        'repetition': equals(),
        'set': equals(),
        'acquisition_time_stamp': ignore(),
        'physiology_time_stamp': each(ignore()),
        'image_type': equals(),
        'image_index': equals(),
        'image_series_index': ignore(),
        'user_int': each(equals()),
        'user_float': each(approx()),
        'attribute_string_len': ignore()
    }

    def check_image_header(output, reference):

        if not output:
            pytest.fail("Missing output")

        if not reference:
            pytest.fail("Missing reference")

        output = output.getHead()
        reference = reference.getHead()

        for attribute, rule in header_rules.items():
            if not rule(getattr(output, attribute), getattr(reference, attribute)):
                print(output)
                print(reference)

                pytest.fail(
                    "Image header '{}' does not match reference. [index {}, series {}]".format(
                        attribute,
                        output.image_index,
                        output.image_series_index
                    )
                )

    try:
        with ismrmrd.File(output_file, 'r') as output_file:
            with ismrmrd.File(reference_file, 'r') as reference_file:
                output_images = output_file[output_group].images or []
                reference_images = reference_file[reference_group].images or []

                for output_image, reference_image in itertools.zip_longest(output_images, reference_images):
                    check_image_header(output_image, reference_image)

    except OSError as e:
        pytest.fail(str(e))

    except RuntimeError as e:
        pytest.fail(str(e))


@pytest.fixture
def validate_stdout(tmp_path: Path) -> Callable:
    def _validate_stdout() -> None:
        files = glob.glob(os.path.join(tmp_path, 'gadgetron_worker*.log.out'))
        files.extend(glob.glob(os.path.join(tmp_path, '*gadgetron.log.out')))

        for file in files:
            if os.stat(file).st_size != 0:
                pytest.fail("stdout is not empty as indicated by {}".format(file))
    
    return _validate_stdout

case_path = './cases'

def load_cases() -> List[Dict[str, str]]:    
    case_list = []
    
    for filename in os.listdir(case_path):
        if not filename.endswith(".yml"):
            continue
        
        with open(os.path.join(case_path, filename), 'r') as file:
            cases = yaml.safe_load(file)

            for case in cases['cases']:
                if not 'mode' in case:
                    case['mode'] = ['server', 'stream', 'distributed']

                case_list.append(case)

    case_list.sort(key=lambda x: (x['name']))

    return case_list

test_cases= []
test_ids = []

def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    get_test_cases(metafunc.config.getoption('--mode'))

def get_test_cases(mode: str) -> None:
    global test_cases, test_ids 
    
    raw_test_cases = load_cases()

    for case in raw_test_cases:
        if mode == 'server':
            send = 'send_to_gadgetron'
            start = 'start_gadgetron_sever'

        elif mode == 'stream':
            send = 'stream_to_gadgetron'
            start = 'start_gadgetron_streaming'

        elif mode == 'distributed':
            send = 'send_to_gadgetron'
            start = 'start_gadgetron_sever_with_additional_nodes'

        else:
            send = 'send_to_gadgetron'
            start = 'start_gadgetron_sever'

            if 'distributed' in case:
                start = 'start_gadgetron_sever_with_additional_nodes'

            if 'stream' in case['reconstruction']:    
                send = 'stream_to_gadgetron'
                start = 'start_gadgetron_streaming'

        test_cases.append((case, start, send))
        test_ids.append(case['name'])


@pytest.fixture
def start_gadgetron(request: pytest.FixtureRequest, external: bool, start_gadgetron_sever: Callable, start_gadgetron_sever_with_additional_nodes: Callable, check_requirements: Callable) -> Callable:
    def _start_gadgetron(fileConfig: Dict[str, str]) -> None:
        local = False        
        if request.param == 'start_gadgetron_sever':
            if 'server' not in fileConfig['mode']:
                pytest.skip("test can't be run in server mode")

            print("starting gadgetron")
            start_gadgetron_sever()

        elif request.param == 'start_gadgetron_sever_with_additional_nodes':
            if 'distributed' not in fileConfig['mode']:
                pytest.skip("test can't be run in distributed mode")

            start_gadgetron_sever_with_additional_nodes(fileConfig)

        elif request.param == 'start_gadgetron_streaming':
            if 'stream' not in fileConfig['mode']:
                pytest.skip("test can't be run in stream mode")

            if external:
                pytest.skip("Stream tests do not run in external mode")
            
            local = True

        else:
            pytest.fail("unknown start_gadgetron type")

        check_requirements(fileConfig['requirements'], fileConfig['tags'], local=local)

    yield _start_gadgetron

@pytest.fixture
def process_data(request: pytest.FixtureRequest, send_to_gadgetron: Callable, stream_to_gadgetron: Callable) -> Callable:
    if request.param == 'send_to_gadgetron':
        return send_to_gadgetron

    if request.param == 'stream_to_gadgetron':
        return stream_to_gadgetron

    pytest.fail("unknown process_data type")

@pytest.mark.parametrize('config, start_gadgetron, process_data', test_cases, ids=test_ids, indirect=['process_data', 'start_gadgetron'])
def test_reconstruction(config: Dict[str, str], start_gadgetron: Callable, process_data: Callable, siemens_to_ismrmrd: Callable, validate_output_data: Callable, validate_dataset_output: Callable, validate_stdout: Callable, start_storage: Callable) -> None:
    storage = start_storage()
    
    start_gadgetron(config)

    if 'dependency' in config:
        dependency_file = siemens_to_ismrmrd(config['dependency'], "dependency")
        _ = process_data(config['dependency'], dependency_file, "dependency")

    reconstruction_file = siemens_to_ismrmrd(config['reconstruction'], "reconstruction")
    output_file = process_data(config['reconstruction'], reconstruction_file, "reconstruction")

    validate_stdout()

    if 'equals' in config['validation']:
        validate_dataset_output(config['validation']['equals'])

    for output_image in config['validation']['images']:   
        validate_output_data(config['validation']['images'][output_image], output_file, output_image)             

    if storage != None:
        storage.kill()
