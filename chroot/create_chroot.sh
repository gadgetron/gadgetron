#!/bin/bash

if [ $(id -u) -ne 0 ]; then
  echo -e "\nPlease start the script as a root or sudo!\n"
  exit 1
else
  if [ $# -ge 4 ]; then

# --ARGUMENTS-- (example)

# CHROOT_GADGETRON_INSTALL_PREFIX:    /usr/local/gadgetron
# CHROOT_GADGETRON_BINARY_DIR:        /home/ubuntu/gadgetron/build
# CHROOT_GIT_SHA1_HASH:               f4d7a9189fd21b07e482d28ecb8b07e589f81f9e
# CHROOT_LIBRARY_PATHS:               /usr/local/lib:/usr/lib/x86_64-linux-gnu
# PACKAGES_PATH:                      /home/ubuntu/packages
    CHROOT_GADGETRON_INSTALL_PREFIX=${1}
    echo CHROOT_GADGETRON_INSTALL_PREFIX: ${CHROOT_GADGETRON_INSTALL_PREFIX}
    CHROOT_GADGETRON_BINARY_DIR=${2}
    echo CHROOT_GADGETRON_BINARY_DIR: ${CHROOT_GADGETRON_BINARY_DIR}
    CHROOT_GIT_SHA1_HASH=${3}
    echo CHROOT_GIT_SHA1_HASH: ${CHROOT_GIT_SHA1_HASH}
    CHROOT_LIBRARY_PATHS=${4}
    echo CHROOT_LIBRARY_PATHS: ${CHROOT_LIBRARY_PATHS}
 
    # Add LIBRARY_PATHS to LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CHROOT_LIBRARY_PATHS}
    export LC_ALL=C
    echo "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}"

    rm -rf ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root
    mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root
    touch ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/source-manifest.txt
    echo "gadgetron    ${CHROOT_GIT_SHA1_HASH}" > ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/source-manifest.txt
    mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron

    if [ $# -eq 5 ]; then
      PACKAGES_PATH=${5}
      echo PACKAGES_PATH: ${PACKAGES_PATH}
      mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/debian/
      cp ${PACKAGES_PATH}/*.deb ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/debian
#      OPENBLAS_PACKAGE=$(basename ${PACKAGES_PATH}/*openblas*.deb)
#      echo "OPENBLAS_PACKAGE: ${OPENBLAS_PACKAGE}"
#      ISMRMRD_PACKAGE=$(basename ${PACKAGES_PATH}/ismrmrd*.deb)
#      echo "ISMRMRD_PACKAGE: ${ISMRMRD_PACKAGE}"
#      CONVERTER_PACKAGE=$(basename ${PACKAGES_PATH}/*siemens*.deb)
#      echo "CONVERTER_PACKAGE: ${CONVERTER_PACKAGE}"
#      GADGETRON_PACKAGE=$(basename ${PACKAGES_PATH}/*gadgetron*.deb)
#      echo "GADGETRON_PACKAGE: ${GADGETRON_PACKAGE}"
    fi

    mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups

    apt-get install debootstrap -y

    debootstrap --variant=buildd --arch amd64 trusty ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron http://gb.archive.ubuntu.com/ubuntu/

    chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get --yes install software-properties-common

    if [ $# -eq 5 ]; then
      cd ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/debian
      dpkg-scanpackages . /dev/null | gzip -9c > Packages.gz
      sed -i '1s;^;deb file:/debian ./\n;' ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/etc/apt/sources.list
    fi

    if [ $# -eq 4 ]; then
      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron add-apt-repository "deb http://gadgetronubuntu.s3.amazonaws.com trusty main"
    fi

    chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron add-apt-repository "deb http://us-east-1.ec2.archive.ubuntu.com/ubuntu/ trusty restricted main multiverse universe"
    chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron add-apt-repository "deb http://us-east-1.ec2.archive.ubuntu.com/ubuntu/ trusty-updates universe restricted multiverse main"
    chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron add-apt-repository "deb http://security.ubuntu.com/ubuntu trusty-security main universe"
    chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron add-apt-repository "deb http://gb.archive.ubuntu.com/ubuntu trusty main"

    chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get update
    chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get --yes install sudo build-essential python-dev python-twisted python-psutil python-numpy python-h5py gdebi-core libxml2-dev

    if [ $# -eq 4 ]; then
      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get --yes --force-yes install libopenblas-base ismrmrd siemens-to-ismrmrd gadgetron
    fi

    if [ $# -eq 5 ]; then
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi $OPENBLAS_PACKAGE
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi $ISMRMRD_PACKAGE
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi $CONVERTER_PACKAGE
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi $GADGETRON_PACKAGE


#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi libopenblas-base_0.2.8-6ubuntu2_amd64.deb
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi libismrmrd1.2_1.2.1-3_amd64.deb
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi ismrmrd-schema_1.2.1-3_all.deb
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi libismrmrd-dev_1.2.1-3_amd64.deb
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi ismrmrd-tools_1.2.1-3_amd64.deb
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi libgadgetron3.1_3.1.1-1_amd64.deb
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi gadgetron-config_3.1.1-1_amd64.deb
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi gadgetron-schema_3.1.1-1_amd64.deb
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi gadgetron-chroot-scripts_3.1.1-1_amd64.deb
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi gadgetron-python-scripts_3.1.1-1_amd64.deb
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi libgadgetron-dev_3.1.1-1_amd64.deb
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi gadgetron-tools_3.1.1-1_amd64.deb
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi gadgetron-all_3.1.1-1_amd64.deb
#      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi siemens-to-ismrmrd_1.0.0-8_amd64.deb
      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi /debian/libopenblas-base_0.2.8-6ubuntu2_amd64.deb
      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi /debian/gadgetron-all_3.1.1-1_amd64.deb
      chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron gdebi /debian/siemens-to-ismrmrd_1.0.0-8_amd64.deb
    fi

    cp -n ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron${CHROOT_GADGETRON_INSTALL_PREFIX}/config/gadgetron.xml.example ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron${CHROOT_GADGETRON_INSTALL_PREFIX}/config/gadgetron.xml

    TAR_FILE_NAME=gadgetron-`date '+%Y%m%d-%H%M'`-${CHROOT_GIT_SHA1_HASH:0:8}
    IMAGE_FILE_NAME=${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/${TAR_FILE_NAME}.img

    tar -zcf "${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/${TAR_FILE_NAME}.tar.gz" --directory "${CHROOT_GADGETRON_BINARY_DIR}/chroot" --exclude=./chroot-root/gadgetron/dev --exclude=./chroot-root/gadgetron/sys --exclude=./chroot-root/gadgetron/proc --exclude=./chroot-root/gadgetron/root ./chroot-root

    dd if=/dev/zero of=${IMAGE_FILE_NAME} bs=3072k seek=1024 count=0
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
    echo -e "\nUsage:  $0 (gadgetron install prefix) (gadgetron binary dir) (GADGETRON_GIT_SHA1_HASH) (LIBRARY_PATHS) optional:(PACKAGES_PATH)\n"
    exit 1
  fi
fi
