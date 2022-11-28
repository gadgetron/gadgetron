#!/bin/bash

set -euo pipefail
export LANG=C

cd "${PREFIX}/share/gadgetron/test/integration" || exit 1

python get_data.py

# confirm reported gadgetron configuration
gadgetron --info

if [[ $(uname) =~ Darwin ]]; then
   echo "Reported hardware and RAM for macOS:"
   system_profiler SPHardwareDataType
   system_profiler SPMemoryDataType
   sysctl hw.memsize
   sw_vers
   # echo "Tests for macOS/Darwin TBD"
   python run_tests.py cases/*
else
   python run_tests.py --ignore-requirements python,cuda cases/*
fi
