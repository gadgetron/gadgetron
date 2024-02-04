#!/bin/bash

set -euo pipefail
export LANG=C

cd "${PREFIX}/share/gadgetron/test/e2e" || exit 1

pytest --ignore-requirements python,cuda
