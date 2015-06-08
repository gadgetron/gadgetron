#!/bin/bash

if [ $(id -u) -ne 0 ]; then
    echo -e "\nPlease start the script as a root or sudo!\n"
    exit 1
else
    if [ $# -ge 1 ]; then

        # --ARGUMENTS-- (example)

        # CHROOT_GADGETRON_BINARY_DIR:        /home/ubuntu/gadgetron/build

        # -----------------------------------------------------------------------------------
        # input parameters
        # -----------------------------------------------------------------------------------

        CHROOT_GADGETRON_BINARY_DIR=${1}
        echo CHROOT_GADGETRON_BINARY_DIR: ${CHROOT_GADGETRON_BINARY_DIR}

        rm -rf ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root
        mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root
        mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups

        # install debian bootstrap
        apt-get install debootstrap -y
        debootstrap --variant=buildd --arch amd64 trusty ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron http://gb.archive.ubuntu.com/ubuntu/

        # install python libraries
        chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get install software-properties-common python-dev python-twisted python-psutil python-numpy python-libxml2 -y
        chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron add-apt-repository "deb http://us.archive.ubuntu.com/ubuntu/ trusty main restricted multiverse universe"  
        chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get update  
        chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get install python-h5py libhdf5-serial-dev hdf5-tools python-pip libplplot-dev -y 

        TAR_FILE_NAME=gadgetron-base-`date '+%Y%m%d-%H%M'`

        tar -zcf "${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/${TAR_FILE_NAME}.tar.gz" --directory "${CHROOT_GADGETRON_BINARY_DIR}/chroot" --exclude=./chroot-root/gadgetron/var --exclude=./chroot-root/gadgetron/dev --exclude=./chroot-root/gadgetron/sys --exclude=./chroot-root/gadgetron/proc --exclude=./chroot-root/gadgetron/root ./chroot-root

        rm -rf "${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root"

        chmod 666 "${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/${TAR_FILE_NAME}.tar.gz"
        exit 0
    else
        echo -e "\nUsage:  $0 (gadgetron binary dir)\n"
        exit 1
    fi
fi
