#!/bin/bash

set -eu

WORKDIR="$(readlink -f $(dirname "$0")/dep-build)"
mkdir -p "$WORKDIR"

PACKAGE_PATH="${WORKDIR}/package"
mkdir -p "$PACKAGE_PATH"

#ranges-v3 (not included with conda)
cd "${WORKDIR}" && \
    rm -rf range-v3 && \
    git clone https://github.com/ericniebler/range-v3.git && \
    cd range-v3 && \
    cp -r include "$PACKAGE_PATH"/

#ISMRMRD
cd "${WORKDIR}" && \
    rm -rf ismrmrd && \
    git clone https://github.com/ismrmrd/ismrmrd.git && \
    cd ismrmrd && \
    mkdir build && \
    cd build && \
    cmake -G Ninja -DCMAKE_INSTALL_PREFIX="$PACKAGE_PATH" ../ && \
    ninja install

#SIEMENS_TO_ISMRMRD
cd "${WORKDIR}" && \
    rm -rf siemens_to_ismrmrd && \
    git clone https://github.com/ismrmrd/siemens_to_ismrmrd.git && \
    cd siemens_to_ismrmrd && \
    mkdir build && \
    cd build && \
    cmake -G Ninja -DBUILD_DYNAMIC=ON -DCMAKE_PREFIX_PATH="$PACKAGE_PATH" -DCMAKE_INSTALL_PREFIX="$PACKAGE_PATH" ../ && \
    ninja install

rsync -a "${PACKAGE_PATH}/" "$CONDA_PREFIX"