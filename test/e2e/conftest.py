#!/usr/bin/env python3

import pytest
import os
import glob
import shutil

from pathlib import Path
from typing import List, Set


def pytest_exception_interact(node, call, report):
    if report.failed and node.config.getoption('--echo-log-on-failure'):
        if hasattr(node, 'funcargs') and 'tmp_path' in node.funcargs:
            tmp_path = node.funcargs['tmp_path']
            for log in glob.glob(os.path.join(tmp_path, '*.log*')):
                with open(log, 'r') as logfile:
                    logdata = logfile.read()
                    report.sections.append((log, logdata))

def pytest_runtest_teardown(item, nextitem):
    if item.config.getoption('--save-results'):
        output_path = item.config.getoption('--save-results')
        output_path = os.path.join(os.path.abspath(output_path), item.callspec.id)
        if hasattr(item, 'funcargs') and 'tmp_path' in item.funcargs:
            tmp_path = item.funcargs['tmp_path']
            shutil.copytree(tmp_path, output_path, dirs_exist_ok=True)


def pytest_addoption(parser):
    parser.addoption(
        '--data-host', action='store', default='https://gadgetronmrd2testdata.blob.core.windows.net/gadgetronmrd2testdata/',
        help='Host from which to download test data.'
    )
    parser.addoption(
        '--cache-disable', action='store_true', default=False,
        help='Disables local caching of input files.'
    )
    parser.addoption(
        '--cache-path', action='store', default=os.path.join(os.path.dirname(__file__), "data"),
        help='Location for storing cached data files.'
    )
    parser.addoption(
        '--ignore-requirements', action='store', nargs='+',
        help="Run tests regardless of whether Pingvin has the specified capabilities."
    )
    parser.addoption(
        '--tags', action='store', nargs='+',
        help='Only run tests with the specified tags.'
    )
    parser.addoption(
        '--echo-log-on-failure', action='store_true', default=False,
        help='Capture and print Pingvin logs on test failure.'
    )
    parser.addoption(
        '--save-results', action='store', default="",
        help='Save Pingvin output and logs in the specified directory.'
    )


@pytest.fixture
def data_host_url(request) -> str:
    return request.config.getoption('--data-host')

@pytest.fixture
def cache_disable(request) -> bool:
    return request.config.getoption('--cache-disable')

@pytest.fixture
def cache_path(request, cache_disable) -> Path:
    if cache_disable:
        return None
    return Path(os.path.abspath(request.config.getoption('--cache-path')))

@pytest.fixture
def ignore_requirements(request) -> Set[str]:
    reqs = request.config.getoption('--ignore-requirements')
    if not reqs:
        return set()
    return set(reqs)

@pytest.fixture
def run_tags(request) -> Set[str]:
    tags = request.config.getoption('--tags')
    if not tags:
        return set()
    return set(tags)
