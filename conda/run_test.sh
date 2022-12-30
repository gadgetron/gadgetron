#!/bin/bash

set -euo pipefail
export LANG=C

cd "${PREFIX}/share/gadgetron/test/integration" || exit 1

python get_data.py

# confirm reported gadgetron configuration
gadgetron --info

if [[ $(uname) =~ Darwin ]]; then
   echo "Place for running unit tests on macOS."
   # test_all
   # echo "Now run integration tests."
   # python run_tests.py cases/*
else
   python run_tests.py --ignore-requirements python,cuda cases/*
fi
