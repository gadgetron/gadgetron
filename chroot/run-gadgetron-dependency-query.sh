#!/bin/bash

BASEDIR=$(dirname $0)

if [ $(id -u) -ne 0 ]; then
  echo -e "\nPlease start the script as a root or sudo!\n"
  exit 1
else
  if [ $# -ge 4 ]; then

    MOUNT_POINT=${1}
    GT_HOST=${2}
    GT_PORT=${3}
    QUERY_OUT=${4}

    if [ $# -eq 5 ]; then
        FULL_PATH_TO_IMG_FILE=${5}

        if find "${MOUNT_POINT}/chroot-root/gadgetron" -maxdepth 0 -empty | read v; then
            mkdir -p ${MOUNT_POINT}
            mount -o loop ${FULL_PATH_TO_IMG_FILE} ${MOUNT_POINT}
        fi
    fi

    chroot ${MOUNT_POINT}/chroot-root/gadgetron /gadgetron-dependency-query.sh $GT_HOST $GT_PORT $QUERY_OUT
    exit $?
  else
    echo -e "\nUsage: $0 <mount point> <Host> <port> <query out file> <optional: full path to img file>\n"
    exit 1
  fi
fi
