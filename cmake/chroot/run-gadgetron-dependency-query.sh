#!/bin/bash

BASEDIR=$(dirname $0)

if [ $(id -u) -ne 0 ]; then
  echo -e "\nPlease start the script as a root or sudo!\n"
  exit 1
else
  if [ $# -ge 6 ]; then

    MOUNT_POINT=${1}
    GT_HOST=${2}
    GT_PORT=${3}
    QUERY_OUT=${4}
    TIME_OUT=${5}
    CHROOT_TIME_OUT=${6}

    if [ $# -eq 7 ]; then
        FULL_PATH_TO_IMG_FILE=${7}

        if find "${MOUNT_POINT}/chroot-root/gadgetron" -maxdepth 0 -empty | read v; then
            mkdir -p ${MOUNT_POINT}
            mount -o loop ${FULL_PATH_TO_IMG_FILE} ${MOUNT_POINT}
        fi
    fi

    if hash timeout 2>/dev/null; then
        ecode=124
        try_count=0
        while [ "$ecode" -eq 124 ] && [ "$try_count" -lt  10 ]; do
            nic=$(timeout $CHROOT_TIME_OUT chroot ${MOUNT_POINT}/chroot-root/gadgetron /home/javeda2/anaconda3/envs/gadgetron/share/gadgetron/chroot/gadgetron-dependency-query.sh $GT_HOST $GT_PORT $QUERY_OUT $TIME_OUT)
            ecode=$?
            try_count=`expr $try_count + 1`
        done
    else
        chroot ${MOUNT_POINT}/chroot-root/gadgetron /home/javeda2/anaconda3/envs/gadgetron/share/gadgetron/chroot/gadgetron-dependency-query.sh $GT_HOST $GT_PORT $QUERY_OUT $TIME_OUT
    fi

    exit $?
  else
    echo -e "\nUsage: $0 <mount point> <Host> <port> <query out file> <connection time out in ms> <script time out in second> <optional: full path to img file>\n"
    exit 1
  fi
fi
