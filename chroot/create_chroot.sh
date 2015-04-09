#!/bin/bash

if [ $(id -u) -ne 0 ]; then
  echo -e "\nPlease start the script as a root or sudo!\n"
  exit 1
else
  if [ $# -ge 5 ]; then

    # --ARGUMENTS-- (example)

    # CHROOT_GADGETRON_INSTALL_PREFIX:    /usr/local/gadgetron
    # CHROOT_GADGETRON_BINARY_DIR:        /home/ubuntu/gadgetron/build
    # CHROOT_LIBRARY_PATHS:               /usr/local/lib:/usr/lib/x86_64-linux-gnu
    # CHROOT_CUDA_LIBRARY:                
    # CHROOT_GADGETRON_SOURCE_DIR:        /home/ubuntu/gadgetron

    CHROOT_GADGETRON_INSTALL_PREFIX=${1}
    echo CHROOT_GADGETRON_INSTALL_PREFIX: ${CHROOT_GADGETRON_INSTALL_PREFIX}
    CHROOT_GADGETRON_BINARY_DIR=${2}
    echo CHROOT_GADGETRON_BINARY_DIR: ${CHROOT_GADGETRON_BINARY_DIR}
    CHROOT_LIBRARY_PATHS=${3}
    echo CHROOT_LIBRARY_PATHS: ${CHROOT_LIBRARY_PATHS}
    CHROOT_CUDA_LIBRARY=${4}
    echo CHROOT_CUDA_LIBRARY: ${CHROOT_CUDA_LIBRARY}
    CHROOT_GADGETRON_SOURCE_DIR=${5}
    echo CHROOT_GADGETRON_SOURCE_DIR: ${CHROOT_GADGETRON_SOURCE_DIR}

    CHROOT_IMAGE_SIZE=1536
    echo CHROOT_IMAGE_SIZE: ${CHROOT_IMAGE_SIZE} MB

    if [ $# -ge 6 ]; then
        CHROOT_SIEMENS_TO_ISMRMRD_EXE=${6} 
        echo CHROOT_SIEMENS_TO_ISMRMRD_EXE: ${CHROOT_SIEMENS_TO_ISMRMRD_EXE}
    else
      echo "SIEMENS_TO_ISMRMRD_EXE not set"
    fi

    # ----------------------------------------------------------------------------------------

    # Add LIBRARY_PATHS to LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CHROOT_LIBRARY_PATHS}
    export LC_ALL=C
    echo "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}"

    # create folders and manifest file
    rm -rf ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root
    mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root
    touch ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/source-manifest.txt

    # get gadgetron SHA1 key
    GADGETRON_INFO=${CHROOT_GADGETRON_INSTALL_PREFIX}/bin/gadgetron_info
    if [ -f ${GADGETRON_INFO} ]; then
      res=$(${GADGETRON_INFO})
      re=".*-- Git SHA1           : ([0-9a-z]+).*"
      if [[ $res =~ $re ]]; then 
        CHROOT_GIT_SHA1_HASH=${BASH_REMATCH[1]}
      fi
    fi

    echo "gadgetron    ${CHROOT_GIT_SHA1_HASH}" > ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/source-manifest.txt
    mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron
    mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups

    # try to find chroot base package
    CHROOT_BASE=`ls -t ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/gadgetron-base*.tar.gz | head -1`

    if [[ -f "${CHROOT_BASE}" ]]; then
        echo "find gadgetron chroot base package : ${CHROOT_BASE}"
    else       
        echo "Cannot find gadgetron chroot base package"
        echo "Creating chroot base package ... "
        ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/create_chroot_base.sh ${CHROOT_GADGETRON_BINARY_DIR}
        CHROOT_BASE=`ls -t ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/gadgetron-base*.tar.gz | head -1`
    fi

    # create chroot from base
    echo "Creating chroot package from  base " ${CHROOT_BASE}
    ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/create_chroot_from_base.sh ${CHROOT_GADGETRON_INSTALL_PREFIX} ${CHROOT_GADGETRON_BINARY_DIR} ${CHROOT_LIBRARY_PATHS} ${CHROOT_CUDA_LIBRARY} ${CHROOT_GADGETRON_SOURCE_DIR} ${CHROOT_BASE} ${CHROOT_IMAGE_SIZE} ${CHROOT_SIEMENS_TO_ISMRMRD_EXE}

    exit 0
  else
    echo -e "\nUsage:  $0 (gadgetron install prefix) (gadgetron binary dir) (LIBRARY_PATHS) (CHROOT_CUDA_LIBRARY) (gadgetron source dir) (SIEMENS_TO_ISMRMRD_EXE)\n"
    exit 1
  fi
fi
