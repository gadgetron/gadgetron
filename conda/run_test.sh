#!/bin/bash

set -euo pipefail
export LANG=C

cd "${PREFIX}/share/gadgetron/test/integration" || exit 1

python get_data.py
python run_tests.py cases/*
