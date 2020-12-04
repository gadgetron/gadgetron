#!/bin/bash

if [ "$EUID" -ne 0 ]
  then echo "Please run as root"
  exit
fi

set -e

GADGETRON_URL="https://github.com/gadgetron/gadgetron"
GADGETRON_BRANCH="master"
WORKDIR="/opt/code"

usage() {
    echo "build-gadgetron.sh            [--source-repo-url ${GADGETRON_URL}]"
    echo "                              [--source-branch ${GADGETRON_BRANCH}]"
    echo "                              [--work-dir ${WORKDIR}]"
}

# Parse command line arguments
while [ "$1" != "" ]; do
    case $1 in
        -r | --source-repo-url )  shift
                                  GADGETRON_URL=$1
                                  ;;
        -b | --source-branch )    shift
                                  GADGETRON_BRANCH=$1
                                  ;;
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

#GADGETRON
cd ${WORKDIR} && \
    rm -rf gadgetron && \
    git clone ${GADGETRON_URL} --branch ${GADGETRON_BRANCH} --single-branch && \
    cd gadgetron && \
    mkdir build && \
    cd build && \
    cmake ../ -G Ninja && \
    ninja && \
    ninja install && \
    ${WORKDIR}/gadgetron/docker/manifest --key .io.gadgetron.gadgetron.sha1 --value `git rev-parse HEAD` && \
    cp ${WORKDIR}/gadgetron/docker/start_supervisor /opt/ && \
    cp ${WORKDIR}/gadgetron/docker/supervisord.conf /opt/

# Install Python interface.
pip3 install ismrmrd multimethod 
pip3 install gadgetron

#HASH for ISMRMRD
cd ${WORKDIR}/ismrmrd && \
    ${WORKDIR}/gadgetron/docker/manifest --key .io.gadgetron.ismrmrd.sha1 --value `git rev-parse HEAD` 

#SIEMENS_TO_ISMRMRD
cd ${WORKDIR} && \
    rm -rf siemens_to_ismrmrd && \
    git clone https://github.com/ismrmrd/siemens_to_ismrmrd.git && \
    cd siemens_to_ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ -G Ninja && \
    ninja && \
    ninja install && \
    ${WORKDIR}/gadgetron/docker/manifest --key .io.gadgetron.siemens_to_ismrmrd.sha1 --value `git rev-parse HEAD` 

#PHILIPS_TO_ISMRMRD
cd ${WORKDIR} && \
    rm -rf philips_to_ismrmrd && \
    git clone https://github.com/ismrmrd/philips_to_ismrmrd.git && \
    cd philips_to_ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ -G Ninja && \
    ninja && \
    ninja install && \
    ${WORKDIR}/gadgetron/docker/manifest --key .io.gadgetron.philips_to_ismrmrd.sha1 --value `git rev-parse HEAD` 