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
   cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=${PREFIX} -DCMAKE_C_COMPILER=${PREFIX}/bin/clang -DCMAKE_CXX_COMPILER=${PREFIX}/bin/clang++ -DPython3_EXECUTABLE=${PREFIX}/bin/python3 -DBUILD_PYTHON_SUPPORT=ON -DUSE_MKL=ON -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="${PREFIX}" "${SRC_DIR}"
else
   cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCUDA_COMPUTE_CAPABILITY=ALL -DUSE_MKL=ON -DUSE_CUDA=ON -DCMAKE_INSTALL_PREFIX="${PREFIX}" "${SRC_DIR}" -DPython3_EXECUTABLE="${PREFIX}/bin/python3"
fi

ninja && ninja install

if [[ $(uname) =~ Darwin ]]; then
   echo "Run storage server tests at end of build, as this testing binary is not installed."
   ./apps/gadgetron/test/server_tests &
fi

TEST_DIR="${PREFIX}/share/gadgetron/test/"
mkdir -p "${TEST_DIR}"
rsync -a --exclude 'data' "${SRC_DIR}/test/integration" "${TEST_DIR}"
