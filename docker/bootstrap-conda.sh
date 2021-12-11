#!/bin/bash

set -eu

WORKDIR="$(dirname "$0")/../dep-build"
mkdir -p "$WORKDIR"

#ranges-v3 (not included with conda)
cd "${WORKDIR}" && \
    rm -rf range-v3 && \
    git clone https://github.com/ericniebler/range-v3.git && \
    cd range-v3 && \
    mkdir -p build && \
    cd build && \
    cmake -G Ninja -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" ../ && \
    ninja && \
    ninja install 

#ISMRMRD
cd "${WORKDIR}" && \
    rm -rf ismrmrd && \
    git clone https://github.com/ismrmrd/ismrmrd.git && \
    cd ismrmrd && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" -G Ninja ../  && \
    ninja && \
    ninja install

#SIEMENS_TO_ISMRMRD
cd "${WORKDIR}" && \
    rm -rf siemens_to_ismrmrd && \
    git clone https://github.com/ismrmrd/siemens_to_ismrmrd.git && \
    cd siemens_to_ismrmrd && \
    mkdir build && \
    cd build && \
    cmake -G Ninja -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" -DBUILD_DYNAMIC=ON ../ && \
    ninja && \
    ninja install 
