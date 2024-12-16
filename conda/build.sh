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

ninja install
