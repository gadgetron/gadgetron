#!/bin/bash
set -euo pipefail

cd conda && python3 validate_versions.py

PKG_DIR="${SRC_DIR}/build_pkg"
mkdir -p "${PKG_DIR}"
cd "${PKG_DIR}" || exit 1

if [[ $(uname) =~ Darwin ]]; then
   echo "Reported hardware and RAM for macOS:"
   uname -a
   system_profiler SPHardwareDataType
   # system_profiler SPMemoryDataType
   sysctl hw.memsize
   sw_vers
   # identify Python causing trouble in cloud environment when building in Python support into macOS build
   which python
   python --version
   # also determine which cmake is being called
   which cmake
   cmake --version
   export MACOSX_DEPLOYMENT_TARGET=10.15
   cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=${CONDA_PREFIX}/bin/clang -DCMAKE_CXX_COMPILER=${CONDA_PREFIX}/bin/clang++ -DBUILD_PYTHON_SUPPORT=OFF -DUSE_MKL=ON -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="${PREFIX}" "${SRC_DIR}"
else
   cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCUDA_COMPUTE_CAPABILITY=ALL -DUSE_MKL=ON -DUSE_CUDA=ON -DCMAKE_INSTALL_PREFIX="${PREFIX}" "${SRC_DIR}"
fi

ninja && ninja install

if [[ $(uname) =~ Darwin ]]; then
   ./apps/gadgetron/test/server_tests &
   ./test/test_all
fi

TEST_DIR="${PREFIX}/share/gadgetron/test/"
mkdir -p "${TEST_DIR}"
rsync -a --exclude 'data' "${SRC_DIR}/test/integration" "${TEST_DIR}"
