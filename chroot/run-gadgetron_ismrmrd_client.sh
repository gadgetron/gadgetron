#!/bin/bash

BASEDIR=$(dirname $0)

if [ $(id -u) -ne 0 ]; then
  echo -e "\nPlease start the script as a root or sudo!\n"
  exit 1
else
  if [ $# -ge 5 ]; then

    MOUNT_POINT=${1}
    ISMRMRD_FILENAME=${2}
    CONDIG_XML=${3}
    GT_HOST=${4}
    GT_PORT=${5}

    if [ $# -eq 5 ]; then
        FULL_PATH_TO_IMG_FILE=${6}

        if find "${MOUNT_POINT}/chroot-root/gadgetron" -maxdepth 0 -empty | read v; then
            mkdir -p ${MOUNT_POINT}
            mount -o loop ${FULL_PATH_TO_IMG_FILE} ${MOUNT_POINT}
        fi
    fi

    chroot ${MOUNT_POINT}/chroot-root/gadgetron /gadgetron_ismrmrd_client.sh $ISMRMRD_FILENAME $CONDIG_XML $GT_HOST $GT_PORT
    exit $?
  else
    echo -e "\nUsage: $0 <mount point> <ismrmrd filename> <config filename> <host> <port> <optional: full path to img file>\n"
    exit 1
  fi
fi
