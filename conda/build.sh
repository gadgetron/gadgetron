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
   export XCODE_ROOT="/Applications/Xcode.app_14.0.1.app"
   # See where else on system python3 may be installed
   # echo "*** Trying to find other Python installations on host system."
   # find /usr -name "python3*" -print
   # find $XCODE_ROOT -name "python3*" -print
   export SYSTEM_DIRECTORIES_TO_IGNORE="${XCODE_ROOT}/Contents/Developer/usr/bin"
   export SYSTEM_DIRECTORIES_TO_IGNORE="${SYSTEM_DIRECTORIES_TO_IGNORE};${XCODE_ROOT}/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX12.3.sdk/usr/lib"
   export SYSTEM_DIRECTORIES_TO_IGNORE="${SYSTEM_DIRECTORIES_TO_IGNORE};${XCODE_ROOT}/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX12.3.sdk/usr/bin"
   export SYSTEM_DIRECTORIES_TO_IGNORE="${SYSTEM_DIRECTORIES_TO_IGNORE};${XCODE_ROOT}/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX12.3.sdk/System/Library/Frameworks"
   export SYSTEM_DIRECTORIES_TO_IGNORE="${SYSTEM_DIRECTORIES_TO_IGNORE};/usr/bin"
   export SYSTEM_DIRECTORIES_TO_IGNORE="${SYSTEM_DIRECTORIES_TO_IGNORE};/usr/lib"
   export SYSTEM_DIRECTORIES_TO_IGNORE="${SYSTEM_DIRECTORIES_TO_IGNORE};/usr/local/bin"
   export SYSTEM_DIRECTORIES_TO_IGNORE="${SYSTEM_DIRECTORIES_TO_IGNORE};/usr/local/lib"
   export SYSTEM_DIRECTORIES_TO_IGNORE="${SYSTEM_DIRECTORIES_TO_IGNORE};/usr/local/opt"
   export SYSTEM_DIRECTORIES_TO_IGNORE="${SYSTEM_DIRECTORIES_TO_IGNORE};/usr/local/miniconda/bin"
   cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_IGNORE_PATH=$SYSTEM_DIRECTORIES_TO_IGNORE -DCMAKE_C_COMPILER=${PREFIX}/bin/clang -DCMAKE_CXX_COMPILER=${PREFIX}/bin/clang++ -DBUILD_PYTHON_SUPPORT=OFF -DUSE_MKL=ON -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="${PREFIX}" "${SRC_DIR}"
else
   cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCUDA_COMPUTE_CAPABILITY=ALL -DUSE_MKL=ON -DUSE_CUDA=ON -DCMAKE_INSTALL_PREFIX="${PREFIX}" "${SRC_DIR}"
fi

ninja && ninja install

if [[ $(uname) =~ Darwin ]]; then
   echo "Run storage server tests at end of build, as this testing binary is not installed."
   ./apps/gadgetron/test/server_tests &
fi

TEST_DIR="${PREFIX}/share/gadgetron/test/"
mkdir -p "${TEST_DIR}"
rsync -a --exclude 'data' "${SRC_DIR}/test/integration" "${TEST_DIR}"
