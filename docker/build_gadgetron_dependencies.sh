#!/bin/bash

if [ "$EUID" -ne 0 ]
  then echo "Please run as root"
  exit
fi

set -e

WORKDIR="/opt/code"

usage() {
    echo "build-gadgetron.sh            [--work-dir ${WORKDIR}]"
}

# Parse command line arguments
while [ "$1" != "" ]; do
    case $1 in
        -w | --work-dir )         shift
                                  WORKDIR=$1
                                  ;;
        -h | --help )             usage
                                  exit
                                  ;;
        * )                       usage
                                  exit 1
    esac
    shift
done

mkdir -p ${WORKDIR}

#ISMRMRD
cd ${WORKDIR} && \
    rm -rf ismrmrd && \
    git clone https://github.com/ismrmrd/ismrmrd.git && \
    cd ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../  -G Ninja && \
    ninja && \
    ninja install

#SIEMENS_TO_ISMRMRD
cd ${WORKDIR} && \
    rm -rf siemens_to_ismrmrd && \
    git clone https://github.com/ismrmrd/siemens_to_ismrmrd.git && \
    cd siemens_to_ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ -G Ninja && \
    ninja && \
    ninja install 

#PHILIPS_TO_ISMRMRD
cd ${WORKDIR} && \
    rm -rf philips_to_ismrmrd && \
    git clone https://github.com/ismrmrd/philips_to_ismrmrd.git && \
    cd philips_to_ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ -G Ninja && \
    ninja && \
    ninja install