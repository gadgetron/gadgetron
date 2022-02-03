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

#cpr (not included anywhere)
cd "${WORKDIR}" && \
    rm -rf cpr && \
    git clone -b 1.7.2 https://github.com/libcpr/cpr.git && \
    cd cpr && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX="$PACKAGE_PATH" -DCPR_FORCE_USE_SYSTEM_CURL=True ../ && \
    make install

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

#mrd-storage-server
# Gadgetron relies heavily on setting time-to-live on items stored. Everything works fine on the current
# mrd-storage-server master branch, but every stored item is kept around forever. Swap this to the master
# branch (or perhaps package a conda package with the mrd-storage-server) when you get the chance.
# Consider this a proof of concept. -klk, Feb 3, 2022
cd "${WORKDIR}" && \
  rm -rf mrd-storage-server && \
  git clone -b respect-ttl https://github.com/kristofferknudsen/mrd-storage-server.git && \
  cd mrd-storage-server && \
  env GOPATH="$PACKAGE_PATH" go install

rsync -a "${PACKAGE_PATH}/" "$CONDA_PREFIX"