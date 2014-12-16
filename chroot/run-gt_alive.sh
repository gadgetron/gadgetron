#!/bin/bash

BASEDIR=$(dirname $0)

if [ $(id -u) -ne 0 ]; then
  echo -e "\nPlease start the script as a root or sudo!\n"
  exit 1
else
  if [ $# -ge 3 ]; then

    MOUNT_POINT=${1}
    HOSTNAME=${2}
    PORT=${3}

    if [ $# -eq 4 ]; then
        FULL_PATH_TO_IMG_FILE=${4}

        if find "${MOUNT_POINT}/chroot-root/gadgetron" -maxdepth 0 -empty | read v; then
            mkdir -p ${MOUNT_POINT}
            mount -o loop ${FULL_PATH_TO_IMG_FILE} ${MOUNT_POINT}
        fi
    fi

    chroot ${MOUNT_POINT}/chroot-root/gadgetron /gt_alive.sh $HOSTNAME $PORT
    exit $?
  else
    echo -e "\nUsage: $0 <mount point> <host> <port> <optional: full path to img file>\n"
    exit 1
  fi
fi
