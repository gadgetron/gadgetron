#!/usr/bin/env bash
set -euo pipefail

mkdir -p build-conda
cd build-conda

cmake -GNinja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCUDA_COMPUTE_CAPABILITY=ALL \
    -DUSE_MKL=ON \
    -DUSE_CUDA=ON \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    ../

ninja install
