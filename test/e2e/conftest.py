#!/usr/bin/python3

import pytest
import os
import glob
import shutil

from pathlib import Path
from typing import List

def pytest_exception_interact(node, call, report):
    if report.failed and node.config.getoption('--echo-log-on-failure'):
        if 'tmp_path' in node.funcargs:
            tmp_path = node.funcargs['tmp_path']
            for log in glob.glob(os.path.join(tmp_path, '*.log*')):
                with open(log, 'r') as logfile:
                    logdata = logfile.read()                
                    report.sections.append((log, logdata))       

def pytest_runtest_teardown(item, nextitem):
    if item.config.getoption('--save-results'):
        output_path=item.config.getoption('--save-results')
        output_path = os.path.join(os.path.abspath(output_path), item.callspec.id)

        if 'tmp_path' in item.funcargs:
            tmp_path = item.funcargs['tmp_path']

            shutil.copytree(tmp_path, output_path, dirs_exist_ok=True)
    

def pytest_addoption(parser):
    parser.addoption(
        '--host', action='store', default='localhost', 
        help='Address of (external) Gadgetron host.'
    )
    parser.addoption(
        '--port', action='store', default='9003', 
        help='Port used by Gadgetron.'
    )
    parser.addoption(
        '--storage-port', action='store', default='9113', 
        help='Port used by Gadgetron Storage Server.'
    )
    parser.addoption(
        '--external', action='store_true', default=False,
        help="External, do not start Gadgetron"
    )
    parser.addoption(
        '--data-host', action='store', default='http://gadgetrondata.blob.core.windows.net/gadgetrontestdata/', 
        help='Host from which to download the data.'
    )
    parser.addoption(
        '--cache-disable', action='store_true', default=False, 
        help='Disables local caching of input files.'
    )
    parser.addoption(
        '--cache-path', action='store', default="", 
        help='Location for storing cached data files.'
    )
    parser.addoption(
        '--ignore-requirements', action='store', default="",
        help="Run tests with the specified tags regardless of Gadgetron capabilities."
    )
    parser.addoption(
        '--tag', action='store', default="", 
        help='Only run tests that has the provided tag.'
    )
    parser.addoption(
        '--echo-log-on-failure', action='store_true', default=False, 
        help='capture test logs on a failed test.'
    )
    parser.addoption(
        '--save-results', action='store', default="",
        help='Save Gadgetron output and client logs to specified folder'
    )
    parser.addoption(
        '--mode', action='store', default="", 
        help='run all tests in a specific mode (server, stream, distributed). Defaults to server unless the test specifies another mode'
    )


@pytest.fixture(scope="module")
def host_url(request) -> str:
    return request.config.getoption('--host')

@pytest.fixture(scope="module")
def port(request) -> int:
    return int(request.config.getoption('--port'))

@pytest.fixture(scope="module")
def storage_port(request) -> int:
    return int(request.config.getoption('--storage-port'))

@pytest.fixture(scope="module")
def external(request) -> bool:
    return request.config.getoption('--external')

@pytest.fixture
def ignore_requirements(request) -> List[str]:
    return request.config.getoption('--ignore-requirements').split(',')

@pytest.fixture
def data_host_url(request) -> str:
    return request.config.getoption('--data-host')

@pytest.fixture
def cache_disable(request) -> bool:
    return request.config.getoption('--cache-disable')

@pytest.fixture
def cache_path(request) -> str:
    base = request.config.getoption('--cache-path')

    if base == "":
        current_dir = Path(os.path.dirname(__file__))
        base = os.path.join(current_dir, "data")
    else:
        base = Path(base)

    return base

@pytest.fixture
def run_tag(request) -> str:
    return request.config.getoption('--tag')


