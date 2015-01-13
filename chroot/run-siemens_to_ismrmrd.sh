#!/bin/bash

BASEDIR=$(dirname $0)

if [ $(id -u) -ne 0 ]; then
  echo -e "\nPlease start the script as a root or sudo!\n"
  exit 1
else
  if [ $# -ge 4 ]; then
    
    MOUNT_POINT=${1}
    DAT_FILENAME=${2}
    ISMRMRD_FILENAME=${3}
    SCAN_NO=${4}

    if [ $# -eq 5 ]; then
        FULL_PATH_TO_IMG_FILE=${5}

        if find "${MOUNT_POINT}/chroot-root/gadgetron" -maxdepth 0 -empty | read v; then
            mkdir -p ${MOUNT_POINT}
            mount -o loop ${FULL_PATH_TO_IMG_FILE} ${MOUNT_POINT}
        fi
    fi

    chroot ${MOUNT_POINT}/chroot-root/gadgetron /siemens_to_ismrmrd.sh $DAT_FILENAME $ISMRMRD_FILENAME $SCAN_NO
    exit $?
  else
    echo -e "\nUsage: $0 <mount point> <dat filename> <ismrmrd filename> <scan number> <optional: full path to img file>\n"
    exit 1
  fi
fi
