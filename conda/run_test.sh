#!/bin/bash

set -euo pipefail
export LANG=C

cd "${PREFIX}/share/gadgetron/test/integration" || exit 1

if [[ $(uname) =~ Darwin ]]; then
   echo "Reported hardware and RAM for macOS:"
   system_profiler SPHardwareDataType
   system_profiler SPMemoryDataType
   sysctl hw.memsize
   echo "Tests for macOS/Darwin TBD"
else
   python get_data.py
   python run_tests.py --ignore-requirements python,cuda cases/*
fi
