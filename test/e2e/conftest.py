#!/usr/bin/python3

import pytest
import os

from pathlib import Path


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
        '--mode', action='store', default="", 
        help='run all tests in a specific mode (server, stream, distributed). Defaults to server unless the test specifies another mode'
    )


@pytest.fixture(scope="module")
def host_url(request):
    return request.config.getoption('--host')

@pytest.fixture(scope="module")
def port(request):
    return request.config.getoption('--port')

@pytest.fixture(scope="module")
def storage_port(request):
    return request.config.getoption('--storage-port')

@pytest.fixture(scope="module")
def external(request):
    return request.config.getoption('--external')

@pytest.fixture
def ignore_requirements(request):
    return request.config.getoption('--ignore-requirements').split(',')

@pytest.fixture
def data_host_url(request):
    return request.config.getoption('--data-host')

@pytest.fixture
def cache_disable(request):
    return request.config.getoption('--cache-disable')

@pytest.fixture
def cache_path(request):
    base = request.config.getoption('--cache-path')

    if base == "":
        current_dir = Path(os.path.dirname(__file__))
        base = current_dir/"data"
    else:
        base = Path(base)

    return base


@pytest.fixture
def run_tag(request):
    return request.config.getoption('--tag')


