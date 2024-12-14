from __future__ import annotations

import re
import os
import glob
import pytest
import hashlib
import itertools
import subprocess
import yaml

import mrd

import numpy
import numpy.typing as npt

import socket
import urllib.error
import urllib.request

from dataclasses import dataclass, field
from typing import Dict, List, Callable, Set, Any


all_test_specs = []
pingvin_capabilities = None


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    """Dynamically generates a test for each test case file"""
    global pingvin_capabilities
    pingvin_capabilities = load_pingvin_capabilities()
    print(f"Pingvin capabilities: {pingvin_capabilities}")
    for spec in load_test_cases():
        all_test_specs.append(spec)


@pytest.mark.parametrize('spec', all_test_specs, ids=lambda s: s.id())
def test_e2e(spec, check_requirements, process_data, validate_output):
    """The main test function for each test case"""
    if spec.dependency is not None:
        process_data(spec.dependency)
    output_file = process_data(spec.reconstruction)
    validate_output(spec.validation, output_file)


@pytest.fixture
def check_requirements(spec: Spec, ignore_requirements: Set[str], run_tags: Set[str]):
    """Checks whether each test case should be run based on Pingvin capabilities and test tags"""

    # Check tags first
    if len(run_tags) > 0 and spec.tags != run_tags:
        pytest.skip("Test missing required tags")
    if 'skip' in spec.tags:
        pytest.skip("Test was marked as skipped")

    # Then check requirements
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
            return lambda value: value is not None and parse_memory(str(target)) <= parse_memory(value)

        def is_positive(value):
            return value is not None and int(value) > 0

        def each(validator):
            return lambda values: all([validator(value) for value in values])

        rules = [
            ('system_memory', lambda req: Rule('memory', has_more_than(req), "Not enough system memory.")),
            ('gpu_support', lambda req: Rule('cuda_support', is_enabled, "CUDA support required.")),
            ('gpu_support', lambda req: Rule('cuda_devices', is_positive, "Not enough CUDA devices.")),
            ('gpu_memory', lambda req: Rule('cuda_memory', each(has_more_than(req)), "Not enough graphics memory."))
        ]

        return [(key, rule(section[key])) for key, rule in rules if key in section]

    rules = rules_from_reqs(spec.requirements)
    for _, rule in rules:
        if rule.capability in ignore_requirements:
            continue
        if not rule.is_satisfied(pingvin_capabilities):
            pytest.skip(rule.message)

@pytest.fixture
def fetch_test_data(cache_path: Path, data_host_url: str, tmp_path: Path) -> Callable:
    """Fetches test data from the remote data host and caches it locally"""
    # If cache_path is disabled, the fetched data will live in the test working directory (tmp_path)
    # PyTest automatically cleans up these directories after 3 runs
    def _fetch_test_data(filename: str, checksum: str) -> str:
        if not cache_path:
            destination = os.path.join(tmp_path, filename)
        else:
            destination = os.path.join(cache_path, filename)

        need_to_fetch = True
        if os.path.exists(destination):
            if not os.path.isfile(destination):
                pytest.fail(f"Destination '{destination}' exists but is not a file")

            if not is_valid(destination, checksum):
                print(f"Destination '{destination}' exists file but checksum does not match... Forcing download")
            else:
                need_to_fetch = False

        if need_to_fetch:
            print(f"Fetching test data: {filename}")
            os.makedirs(os.path.dirname(destination), exist_ok=True)
            url = f"{data_host_url}{filename}"
            urlretrieve(url, destination)

        if not is_valid(destination, checksum):
            pytest.fail(f"Downloaded file '{destination}' does not match checksum")

        return destination

    return _fetch_test_data


@pytest.fixture
def process_data(fetch_test_data, tmp_path):
    """Runs the Pingvin on the input test data, producing an output file."""
    def _process_data(job):
        input_file = fetch_test_data(job.datafile, job.checksum)
        output_file = os.path.join(tmp_path, job.name + ".output.mrd")

        invocations = []
        for args in job.args:
            invocations.append(f"pingvin {args}")
        invocations[0] += f" --input {input_file}"
        invocations[-1] += f" --output {output_file}"

        command = " | ".join(invocations)
        command = ['bash', '-c', command]

        print(f"Run: {command}")

        log_stdout_filename = os.path.join(tmp_path, f"pingvin_{job.name}.log.out")
        log_stderr_filename = os.path.join(tmp_path, f"pingvin_{job.name}.log.err")
        with open(log_stdout_filename, 'w') as log_stdout:
            with open(log_stderr_filename, 'w') as log_stderr:
                result = subprocess.run(command, stdout=log_stdout, stderr=log_stderr, cwd=tmp_path)
                if result.returncode != 0:
                    pytest.fail(f"Pingvin failed with return code {result.returncode}. See {log_stderr_filename} for details.")

        return output_file

    return _process_data

@pytest.fixture
def validate_output(fetch_test_data):
    """Validates each image (data and header) in the output file against the reference file."""
    def _validate_output(spec: Spec.Validation, output_file: str) -> None:
        reference_file = fetch_test_data(spec.reference, spec.checksum)

        reference_images = extract_image_data(reference_file)
        output_images = extract_image_data(output_file)

        for test in spec.image_series_tests:
            ref_data = reference_images[test.image_series]['data']
            out_data = output_images[test.image_series]['data']
            validate_image_data(out_data, ref_data, test.scale_comparison_threshold, test.value_comparison_threshold)

            ref_headers = reference_images[test.image_series]['headers']
            out_headers = output_images[test.image_series]['headers']
            for ref, out in itertools.zip_longest(ref_headers, out_headers):
                validate_image_header(out, ref)

    return _validate_output


def load_pingvin_capabilities() -> Dict[str, str]:
    command = ["pingvin", "--info"]
    res = subprocess.run(command, capture_output=True, text=True)
    if res.returncode != 0:
        pytest.fail(f"Failed to query Pingvin capabilities... {res.args} return {res.returncode}")

    pingvin_info = res.stderr

    value_pattern = r"(?:\s*):(?:\s+)(?P<value>.*)?"

    capability_markers = {
        'version': "Version",
        'build': "Git SHA1",
        'memory': "System Memory size",
        'cuda_support': "CUDA Support",
        'cuda_devices': "CUDA Device count"
    }

    plural_capability_markers = {
        'cuda_memory': "CUDA Device Memory size",
    }

    def find_value(marker):
        pattern = re.compile(marker + value_pattern, re.IGNORECASE)
        match = pattern.search(pingvin_info)
        if match:
            return match['value']
        else:
            return None

    def find_plural_values(marker):
        pattern = re.compile(marker + value_pattern, re.IGNORECASE)
        return [match['value'] for match in pattern.finditer(pingvin_info)]

    capabilities = {key: find_value(marker) for key, marker in capability_markers.items()}
    capabilities.update({key: find_plural_values(marker) for key, marker in plural_capability_markers.items()})

    return capabilities


def load_test_cases() -> List[Dict[str, str]]:
    specs = []
    for filename in glob.glob('cases/*.yml'):
        spec = Spec.fromfile(filename)
        specs.append(spec)
    return sorted(specs, key=lambda s: s.id())


def checksum(file: Path) -> str:
    md5 = hashlib.new('md5')
    with open(file, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            md5.update(chunk)
    return md5.hexdigest()

def is_valid(file: Path, digest: str) -> bool:
    if not os.path.isfile(file):
        return False
    return digest == checksum(file)

def urlretrieve(url: str, filename: str, retries: int = 5) -> str:
    if retries <= 0:
        pytest.fail("Download from {} failed".format(url))
    try:
        with urllib.request.urlopen(url, timeout=60) as connection:
            with open(filename,'wb') as f:
                for chunk in iter(lambda : connection.read(1024*1024), b''):
                    f.write(chunk)
            return connection.headers["Content-MD5"]
    except (urllib.error.URLError, ConnectionResetError, socket.timeout) as exc:
        print("Retrying connection for file {}, reason: {}".format(filename, str(exc)))
        return urlretrieve(url, filename, retries=retries-1)

def extract_image_data(filename: Path) -> Dict[int, Dict[str, Any]]:
    dataset = {}
    with mrd.BinaryMrdReader(filename) as reader:
        _ = reader.read_header()
        for item in reader.read_data():
            if type(item.value).__name__ != "Image":
                continue
            head = item.value.head
            if head.image_series_index not in dataset:
                dataset[head.image_series_index] = []
            dataset[head.image_series_index].append(item.value)
    res = {}
    for series, images in dataset.items():
        res[series] = {}
        res[series]['data'] = numpy.stack([img.data for img in images])
        res[series]['headers'] = map(lambda img: img.head, images)
        res[series]['meta'] = map(lambda img: img.meta, images)
    return res

def validate_image_data(output_data: npt.ArrayLike, reference_data: npt.ArrayLike,
                scale_threshold: float, value_threshold: float) -> None:
    assert output_data.shape == reference_data.shape
    assert output_data.dtype == reference_data.dtype

    ref_data = reference_data.flatten().astype('float32')
    out_data = output_data.flatten().astype('float32')

    norm_diff = numpy.linalg.norm(out_data - ref_data) / numpy.linalg.norm(ref_data)
    scale = numpy.dot(out_data, out_data) / numpy.dot(out_data, ref_data)

    if value_threshold < norm_diff:
        pytest.fail("Comparing values, norm diff: {} (threshold: {})".format(norm_diff, value_threshold))

    if scale_threshold < abs(1 - scale):
        pytest.fail("Comparing image scales, ratio: {} ({}) (threshold: {})".format(scale, abs(1 - scale),
                                                                                        scale_threshold))

def validate_image_header(output_header: mrd.ImageHeader, reference_header: mrd.ImageHeader) -> None:
    def equals():
        # Account for *converted* reference header field having value 0 instead of None
        return lambda out, ref: out == ref or (out is None and ref == 0)

    def approx(threshold=1e-6):
        return lambda out, ref: abs(out - ref) <= threshold

    def ignore():
        return lambda out, ref: True

    def each(rule):
        return lambda out, ref: all(rule(out, ref) for out, ref in itertools.zip_longest(out, ref))

    header_rules = {
        'flags': equals(),
        'measurement_uid': equals(),
        'field_of_view': each(approx()),
        'position': each(approx()),
        'col_dir': each(approx()),
        'line_dir': each(approx()),
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
        'image_series_index': equals(),
        # Ignore user values since ISMRMRD->MRD converted files always have user_int/user_float==[0,0,0,0,0,0,0,0]
        'user_int': each(ignore()),
        'user_float': each(ignore())
    }

    for attribute, rule in header_rules.items():
        if not rule(getattr(output_header, attribute), getattr(reference_header, attribute)):
            pytest.fail(f"Image header '{attribute}' does not match reference"
                f" (series {output_header.image_series_index}, index {output_header.image_index})"
                f" [{getattr(output_header, attribute)} != {getattr(reference_header, attribute)}]")


@dataclass
class Spec():
    """Defines a test case specification"""

    @dataclass
    class Job():
        """Defines a job to be run by Pingvin"""
        name: str
        datafile: str
        checksum: str
        args: List[str]

        @staticmethod
        def fromdict(config: Dict[str, str], name: str) -> Spec.Job:
            if not config:
                return None

            datafile = config['data']
            if not datafile:
                raise ValueError(f"Missing 'data' key in {name} configuration")

            checksum = config['checksum']
            if not checksum:
                raise ValueError(f"Missing 'checksum' key in {name} configuration")

            args = []
            if 'run' in config:
                for run in config['run']:
                    args.append(run['args'])
            else:
                args.append(config['args'])

            return Spec.Job(name=name, datafile=datafile, checksum=checksum, args=args)

    @dataclass
    class ImageSeriesTest():
        """Defines a test for an image series comparison"""
        image_series: int
        scale_comparison_threshold: float
        value_comparison_threshold: float

    @dataclass
    class Validation():
        """Defines a validation test for the output of a job"""
        reference: str
        checksum: str
        image_series_tests: List[ImageSeriesTest]

        @staticmethod
        def fromdict(config: Dict[str, str]) -> Spec.Validation:
            if not config:
                return None
            reference = config['reference']
            if not reference:
                raise ValueError("Missing 'reference' key in validation configuration")
            checksum = config['checksum']
            if not checksum:
                raise ValueError("Missing 'checksum' key in validation configuration")
            tests = config['tests']
            if not tests:
                raise ValueError("Missing 'tests' key in validation configuration")
            if not isinstance(tests, list):
                raise ValueError("Key 'tests' should be a list in validation configuration")

            image_series_tests = []
            for test in tests:
                num = test['image_series']
                st = test.get('scale_comparison_threshold', 0.01)
                vt = test.get('value_comparison_threshold', 0.01)
                image_series_tests.append(
                    Spec.ImageSeriesTest(image_series=num, scale_comparison_threshold=st, value_comparison_threshold=vt)
                )

            return Spec.Validation(reference=reference, checksum=checksum,
                    image_series_tests=image_series_tests)

    name: str
    tags: Set[str] = field(default_factory=set)
    requirements: Dict[str, str] = field(default_factory=dict)

    dependency: Spec.Job = None
    reconstruction: Spec.Job = None
    validation: Spec.Validation = None

    def id(self):
        return f"{self.name}"

    @staticmethod
    def fromfile(filename: str) -> Spec:
        with open(filename, 'r') as file:
            parsed = yaml.safe_load(file)
            name = os.path.relpath(filename)
            spec = Spec(name=name)

            tags = parsed.get('tags', None)
            if not tags:
                tags = []
            if not isinstance(tags, list):
                tags = [tags]
            spec.tags = set(tags)

            requirements = parsed.get('requirements', None)
            if not requirements:
                requirements = {}
            if not isinstance(requirements, dict):
                raise ValueError(f"Invalid requirements in {filename}")
            spec.requirements = requirements

            spec.dependency = Spec.Job.fromdict(parsed.get('dependency', None), 'dependency')
            spec.reconstruction = Spec.Job.fromdict(parsed['reconstruction'], 'reconstruction')
            spec.validation = Spec.Validation.fromdict(parsed.get('validation', None))

            return spec