#!/bin/bash

BASEDIR=$(dirname $0)

if [ $(id -u) -ne 0 ]; then
  echo -e "\nPlease start the script as a root or sudo!\n"
  exit 1
else
  if [ $# -ge 7 ]; then

    MOUNT_POINT=${1}
    OUTPUT_FILENAME=${2}
    CONDIG_XML=${3}
    GT_HOST=${4}
    GT_PORT=${5}
    GT_TIMEOUT=${6}
    CHROOT_TIMEOUT=${7}

    if [ $# -eq 8 ]; then
        FULL_PATH_TO_IMG_FILE=${8}

        if find "${MOUNT_POINT}/chroot-root/gadgetron" -maxdepth 0 -empty | read v; then
            mkdir -p ${MOUNT_POINT}
            mount -o loop ${FULL_PATH_TO_IMG_FILE} ${MOUNT_POINT}
        fi
    fi

    if hash timeout 2>/dev/null; then
        ecode=124
        try_count=0
        while [ "$ecode" -eq 124 ] && [ "$try_count" -lt  5 ]; do
            nic=$(timeout $CHROOT_TIMEOUT chroot ${MOUNT_POINT}/chroot-root/gadgetron /home/javeda2/anaconda3/envs/gadgetron/share/gadgetron/chroot/gadgetron_ismrmrd_client_noise_summary.sh $OUTPUT_FILENAME $CONDIG_XML $GT_HOST $GT_PORT $GT_TIMEOUT)
            ecode=$?
            try_count=`expr $try_count + 1`
        done
    else
        chroot ${MOUNT_POINT}/chroot-root/gadgetron /home/javeda2/anaconda3/envs/gadgetron/share/gadgetron/chroot/gadgetron_ismrmrd_client_noise_summary.sh $OUTPUT_FILENAME $CONDIG_XML $GT_HOST $GT_PORT $GT_TIMEOUT
    fi

    exit $?
  else
    echo -e "\nUsage: $0 <mount point> <output filename> <local config filename> <host> <port> <connection time out in ms> <script time out in second> <optional: full path to img file>\n"
    exit 1
  fi
fi
