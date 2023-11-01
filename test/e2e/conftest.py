#!/usr/bin/python3

import pytest
import os

from pathlib import Path

def pytest_addoption(parser):
    parser.addoption(
        '--host', action='store', default='http://gadgetrondata.blob.core.windows.net/gadgetrontestdata/', 
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


@pytest.fixture
def host_url(request):
    return request.config.getoption('--host')

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
