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

    if [ $# -ge 6 ]; then
        CHROOT_IMAGE_SIZE=${6}         
    else
        CHROOT_IMAGE_SIZE=1536
    fi
    echo CHROOT_IMAGE_SIZE: ${CHROOT_IMAGE_SIZE} MB

    if [ $# -ge 7 ]; then
        CHROOT_SIEMENS_TO_ISMRMRD_EXE=${7} 
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
        tar -xzf ${CHROOT_BASE} -C ${CHROOT_GADGETRON_BINARY_DIR}/chroot
        sleep 3
    else       
        echo "Cannot find gadgetron chroot base package"

        # install debian boot strap
        apt-get install debootstrap -y

        debootstrap --variant=buildd --arch amd64 trusty ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron http://gb.archive.ubuntu.com/ubuntu/

        # install python packages
        chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get install software-properties-common python-dev python-twisted python-psutil python-numpy python-libxml2 -y
        chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron add-apt-repository "deb http://us.archive.ubuntu.com/ubuntu/ trusty main restricted multiverse universe"  
        chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get update  
        chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get install python-h5py libhdf5-serial-dev hdf5-tools python-pip -y 
    fi

    # install gadgetron
    cd ${CHROOT_GADGETRON_BINARY_DIR}
    make install DESTDIR="${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron" -j8

    # This copies the SIEMENS_TO_ISMRMRD executable if it is installed
    if [ $# -ge 7 ]; then      
      cp ${CHROOT_SIEMENS_TO_ISMRMRD_EXE} "${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/${CHROOT_GADGETRON_INSTALL_PREFIX}/bin/"
    fi

    ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/generate_gadgetron_root ${CHROOT_GADGETRON_INSTALL_PREFIX} ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron

    cp ${CHROOT_CUDA_LIBRARY} ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/${CHROOT_GADGETRON_INSTALL_PREFIX}/lib  
    cp -n ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron${CHROOT_GADGETRON_INSTALL_PREFIX}/share/gadgetron/config/gadgetron.xml.example ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron${CHROOT_GADGETRON_INSTALL_PREFIX}/share/gadgetron/config/gadgetron.xml

    # install ismrmrd python
    ISMRMRD_PYTHON_FOLDER=${CHROOT_GADGETRON_INSTALL_PREFIX}/share/gadgetron/python/ismrmrd
    if [ -d ${ISMRMRD_PYTHON_FOLDER} ]; then
      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron pip install cython h5py pyxb
      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron pip install --upgrade h5py
      cp -rf ${ISMRMRD_PYTHON_FOLDER} "${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/${CHROOT_GADGETRON_INSTALL_PREFIX}/share/gadgetron/python"
    fi

    # create tar package
    TAR_FILE_NAME=gadgetron-`date '+%Y%m%d-%H%M'`-${CHROOT_GIT_SHA1_HASH:0:8}
    IMAGE_FILE_NAME=${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/${TAR_FILE_NAME}.img

    tar -zcf "${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/${TAR_FILE_NAME}.tar.gz" --directory "${CHROOT_GADGETRON_BINARY_DIR}/chroot" --exclude=./chroot-root/gadgetron/var --exclude=./chroot-root/gadgetron/dev --exclude=./chroot-root/gadgetron/sys --exclude=./chroot-root/gadgetron/proc --exclude=./chroot-root/gadgetron/root ./chroot-root

    # create image file
    dd if=/dev/zero of=${IMAGE_FILE_NAME} bs=${CHROOT_IMAGE_SIZE}k seek=1024 count=0
    mke2fs -F -t ext3 ${IMAGE_FILE_NAME}
    mkdir ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_root
    mount -o loop ${IMAGE_FILE_NAME} ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_root
    tar -xzf ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/${TAR_FILE_NAME}.tar.gz -C ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_root/
    sleep 3
    umount ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_root
    rmdir ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_root
    rm -rf "${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root"

    chmod 666 "${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/${TAR_FILE_NAME}.tar.gz"
    chmod 666 "${IMAGE_FILE_NAME}"
    exit 0
  else
    echo -e "\nUsage:  $0 (gadgetron install prefix) (gadgetron binary dir) (LIBRARY_PATHS) (CHROOT_CUDA_LIBRARY) (gadgetron source dir) (SIEMENS_TO_ISMRMRD_EXE)\n"
    exit 1
  fi
fi
