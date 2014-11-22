#!/bin/bash

if [ $(id -u) -ne 0 ]; then 
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1

else
 BASEDIR=$(dirname $0)

 if [ $# -eq 2 ]; then

  service gadgetron_chroot stop

  CHROOT_IMAGE_FILENAME=${1}
  echo CHROOT_IMAGE_FILENAME=${CHROOT_IMAGE_FILENAME}

  CHROOT_INSTALL_PATH=${2}
  echo CHROOT_INSTALL_PATH=${CHROOT_INSTALL_PATH}

  mkdir -p ${CHROOT_INSTALL_PATH}

  cp -rf ${CHROOT_IMAGE_FILENAME} ${CHROOT_INSTALL_PATH}/

  FILENAME_WITH_EXTENSION=${CHROOT_IMAGE_FILENAME##*/}
  FILENAME=${FILENAME_WITH_EXTENSION%.*}
  FILENAME=${FILENAME%.*}
  echo ${FILENAME}

  mkdir ${CHROOT_INSTALL_PATH}/${FILENAME}

  echo untar ${CHROOT_INSTALL_PATH}/${FILENAME_WITH_EXTENSION} ... 

  tar -xzf ${CHROOT_INSTALL_PATH}/${FILENAME_WITH_EXTENSION} --directory="${CHROOT_INSTALL_PATH}/${FILENAME}" .

  rm -f ${CHROOT_INSTALL_PATH}/current

  ln -s ${CHROOT_INSTALL_PATH}/${FILENAME} ${CHROOT_INSTALL_PATH}/current

  cp -f ${CHROOT_INSTALL_PATH}/current/chroot-root/gadgetron/webapp/gadgetron_chroot.conf /etc/init/

  service gadgetron_chroot start

  exit 0
 
 else
  echo -e "\nUsage: $0 chroot_image_name chroot_install_path\n"
  exit 1 
 fi
fi
