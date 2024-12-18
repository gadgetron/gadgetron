#!/usr/bin/env bash
set -euo pipefail

mkdir -p build-conda
cd build-conda

if [[ $(uname) =~ Darwin ]]; then
   echo "Reported hardware and RAM for macOS:"
   uname -a
   system_profiler SPHardwareDataType
   # system_profiler SPMemoryDataType
   sysctl hw.memsize
   sw_vers
   cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=${PREFIX} -DCMAKE_C_COMPILER=${CONDA_PREFIX}/bin/clang -DCMAKE_CXX_COMPILER=${CONDA_PREFIX}/bin/clang++ -DPython3_EXECUTABLE=${CONDA_PREFIX}/bin/python3 -DBUILD_PYTHON_SUPPORT=ON -DUSE_MKL=ON -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="${PREFIX}" ..
else
   cmake -GNinja \
       -DCMAKE_BUILD_TYPE=Release \
       -DCUDA_COMPUTE_CAPABILITY=ALL \
       -DUSE_MKL=ON \
       -DUSE_CUDA=ON \
       -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
       -DPython3_EXECUTABLE="${CONDA_PREFIX}/bin/python3" \
       ..
fi

ninja install
