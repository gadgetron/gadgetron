#!/bin/bash
set -euo pipefail

cd conda && python3 validate_versions.py

PKG_DIR="${SRC_DIR}/build_pkg"
mkdir -p "${PKG_DIR}"
cd "${PKG_DIR}" || exit 1

if [[ $(uname) =~ Darwin ]]; then
    cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=${CONDA_PREFIX}/bin/clang -DCMAKE_CXX_COMPILER=${CONDA_PREFIX}/bin/clang++ -DUSE_MKL=ON -DCMAKE_INSTALL_PREFIX="${PREFIX}" "${SRC_DIR}"
else
    cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCUDA_COMPUTE_CAPABILITY=ALL -DUSE_MKL=ON -DUSE_CUDA=ON -DCMAKE_INSTALL_PREFIX="${PREFIX}" "${SRC_DIR}"
fi
ninja && ninja install

TEST_DIR="${PREFIX}/share/gadgetron/test/"
mkdir -p "${TEST_DIR}"
rsync -a --exclude 'data' "${SRC_DIR}/test/integration" "${TEST_DIR}"
