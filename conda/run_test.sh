#!/usr/bin/env bash
set -euo pipefail

cd test/e2e/

test -d cases
test -f conftest.py

pytest