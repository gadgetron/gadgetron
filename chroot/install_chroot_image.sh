#!/bin/bash

function umount_check {
    MAX_TRY=500
    MOUNT_DIR=$1
    UMOUNT_READY=0
    UMOUNT_TRY=0
    while [ ${UMOUNT_READY} -eq 0 ]; do
        if mountpoint -q ${MOUNT_DIR}; then
            let UMOUNT_TRY++
            if [ $UMOUNT_TRY -eq $MAX_TRY ]; then
                UMOUNT_READY=1
            else
                sleep 0.5
            fi
        else
            UMOUNT_READY=1
        fi
    done

    if mountpoint -q ${MOUNT_DIR}; then
        umount ${MOUNT_DIR}
    fi
}

if [ $(id -u) -ne 0 ]; then 
    echo -e "\nPlease start the script as a root or sudo!\n"
    exit 1

else

    BASEDIR=$(dirname $0)

    install_img=0

    if [ $# -eq 1 ]; then
        CHROOT_IMAGE_FILENAME=${1}
        CHROOT_ISMRMRD_DATA_PATH=/tmp/gadgetron_chroot_data
        CHROOT_INSTALL_PATH=/home/gadgetron_chroot
    else
        if [ $# -eq 2 ]; then
            CHROOT_IMAGE_FILENAME=${1}
            CHROOT_ISMRMRD_DATA_PATH=${2}
            CHROOT_INSTALL_PATH=/home/gadgetron_chroot
        else
            if [ $# -eq 3 ]; then
                if [ ${2} == "latest" ]; then
                    TAR_NAME=`find ${1} -type f -name 'gadgetron-*.tar.gz' |sort |head -n1`
                    CHROOT_IMAGE_FILENAME=${TAR_NAME}
                else
                    CHROOT_IMAGE_FILENAME=${1}/${2}
                fi
                CHROOT_INSTALL_PATH=${3}
                CHROOT_ISMRMRD_DATA_PATH=/tmp/gadgetron_chroot_data
            else
                if [ $# -eq 4 ]; then
                    if [ ${2} == "latest" ]; then
                        TAR_NAME=`find ${1} -type f -name 'gadgetron-*.tar.gz' |sort |head -n1`
                        CHROOT_IMAGE_FILENAME=${1}/${TAR_NAME}
                    else
                        CHROOT_IMAGE_FILENAME=${1}/${2}         
                    fi
                    CHROOT_INSTALL_PATH=${3}
                    CHROOT_ISMRMRD_DATA_PATH=${4}
                else
                    if [ $# -eq 5 ]; then
                        if [ ${2} == "latest" ]; then
                            TAR_NAME=`find ${1} -type f -name 'gadgetron-*.tar.gz' |sort |head -n1`
                            CHROOT_IMAGE_FILENAME=${1}/${TAR_NAME}
                        else
                            CHROOT_IMAGE_FILENAME=${1}/${2}         
                        fi
                        CHROOT_INSTALL_PATH=${3}
                        CHROOT_ISMRMRD_DATA_PATH=${4}

                        if [ ${5} -eq 1 ]; then
                            install_img=1
                            IMG_NAME=`find ${1} -type f -name 'gadgetron-*.img' |sort |head -n1`
                            CHROOT_IMAGE_IMG_FILENAME=${IMG_NAME}
                        fi
                    else
                        echo -e "\nUsage 1, install chroot image to /home/gadgetron_chroot: $0 chroot_image_file"
                        echo -e "\nUsage 2, install chroot image to selected install path: $0 chroot_image_file chroot_ismrmrd_data_path"
                        echo -e "\nUsage 3, install chroot image to selected install path: $0 chroot_image_path chroot_image_name chroot_install_path"
                        echo -e "\nUsage 4, : $0 chroot_image_path chroot_image_name chroot_install_path chroot_ismrmrd_data_path"
                        echo -e "\n           install chroot image to selected install path, if chroot_image_name=latest, the newest chroot image in the folder will be installed"
                        echo -e "\nUsage 5, : $0 chroot_image_path chroot_image_name chroot_install_path chroot_ismrmrd_data_path install_img"
                        echo -e "\n           like Usage 4, if install_img=1, the corresponding .img package will be copied to chroot_install_path"
                        exit 1
                    fi
                fi
            fi  
        fi  
    fi

    CHROOT_ISMRMRD_DATA_MOUNT_PATH=/tmp/gadgetron_data

    service gadgetron_chroot stop

    echo CHROOT_IMAGE_FILENAME=${CHROOT_IMAGE_FILENAME}
    echo CHROOT_INSTALL_PATH=${CHROOT_INSTALL_PATH}

    # umount existing folder
    if mountpoint -q ${CHROOT_INSTALL_PATH}/current/chroot-root/gadgetron/proc; then
        umount ${CHROOT_INSTALL_PATH}/current/chroot-root/gadgetron/proc
        umount_check ${CHROOT_INSTALL_PATH}/current/chroot-root/gadgetron/proc
    fi

    if mountpoint -q ${CHROOT_INSTALL_PATH}/current/chroot-root/gadgetron/dev; then
        umount ${CHROOT_INSTALL_PATH}/current/chroot-root/gadgetron/dev
        umount_check ${CHROOT_INSTALL_PATH}/current/chroot-root/gadgetron/dev
    fi

    if mountpoint -q ${CHROOT_INSTALL_PATH}/current/chroot-root/gadgetron/sys; then
        umount ${CHROOT_INSTALL_PATH}/current/chroot-root/gadgetron/sys
        umount_check ${CHROOT_INSTALL_PATH}/current/chroot-root/gadgetron/sys
    fi

    if mountpoint -q ${CHROOT_INSTALL_PATH}/current/chroot-root/gadgetron/tmp/gadgetron_data; then
        umount ${CHROOT_INSTALL_PATH}/current/chroot-root/gadgetron/tmp/gadgetron_data
        umount_check ${CHROOT_INSTALL_PATH}/current/chroot-root/gadgetron/tmp/gadgetron_data
    fi

    # untar the chroot package
    mkdir -p ${CHROOT_INSTALL_PATH}

    cp -rf ${CHROOT_IMAGE_FILENAME} ${CHROOT_INSTALL_PATH}/

    FILENAME_WITH_EXTENSION=${CHROOT_IMAGE_FILENAME##*/}
    FILENAME=${FILENAME_WITH_EXTENSION%.*}
    FILENAME=${FILENAME%.*}
    echo ${FILENAME}

    mkdir ${CHROOT_INSTALL_PATH}/${FILENAME}

    echo untar ${CHROOT_INSTALL_PATH}/${FILENAME_WITH_EXTENSION} ... 
    tar -xzf ${CHROOT_INSTALL_PATH}/${FILENAME_WITH_EXTENSION} --directory="${CHROOT_INSTALL_PATH}/${FILENAME}"

    # detect whether new or older style chroot image is used
    if [ -d ${CHROOT_INSTALL_PATH}/${FILENAME}/chroot-root ]; then
        echo Chroot script generated this package ...
    else
        echo Docker generated this package ...
        rm -rf ${CHROOT_INSTALL_PATH}/${FILENAME}/*
        mkdir ${CHROOT_INSTALL_PATH}/${FILENAME}/chroot-root
        mkdir ${CHROOT_INSTALL_PATH}/${FILENAME}/chroot-root/gadgetron
        tar -xzf ${CHROOT_INSTALL_PATH}/${FILENAME_WITH_EXTENSION} --directory="${CHROOT_INSTALL_PATH}/${FILENAME}/chroot-root/gadgetron"
    fi

    rm -f ${CHROOT_INSTALL_PATH}/current

    ln -s ${CHROOT_INSTALL_PATH}/${FILENAME} ${CHROOT_INSTALL_PATH}/current

    cp -f ${CHROOT_INSTALL_PATH}/current/chroot-root/gadgetron/usr/local/share/gadgetron/chroot/gadgetron_chroot.conf /etc/init/

    # create chroot ismrmrd data mount point outside chroot
    echo CHROOT_ISMRMD_DATA_MOUNT_PATH=$CHROOT_ISMRMRD_DATA_MOUNT_PATH
    echo CHROOT_ISMRMRD_DATA_PATH=$CHROOT_ISMRMRD_DATA_PATH
    rm -rf ${CHROOT_ISMRMRD_DATA_MOUNT_PATH}
    mkdir -p ${CHROOT_ISMRMRD_DATA_PATH}
    cp -f ${CHROOT_INSTALL_PATH}/current/chroot-root/gadgetron/usr/local/share/gadgetron/chroot/clean_gadgetron_data.sh ${CHROOT_ISMRMRD_DATA_PATH}
    ln -s ${CHROOT_ISMRMRD_DATA_PATH} ${CHROOT_ISMRMRD_DATA_MOUNT_PATH}

    # create chroot ismrmrd mount point inside chroot
    if [ ${install_img} -eq 1 ]; then
        echo "copy image file : ${CHROOT_IMAGE_IMG_FILENAME} ... "      
        cp -f ${CHROOT_IMAGE_IMG_FILENAME} ${CHROOT_INSTALL_PATH}/
    fi

    service gadgetron_chroot start

    exit 0
fi
