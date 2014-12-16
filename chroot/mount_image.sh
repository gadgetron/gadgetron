#!/bin/bash

function mount_safe {
  MOUNT_POINT=$1
  MOUNT_DIR=$2
  mkdir -p $MOUNT_POINT
  if find $MOUNT_POINT -maxdepth 0 -empty | read v; then
    mount --bind $MOUNT_DIR $MOUNT_POINT
    MOUNT_READY=0
    MOUNT_TRY=0
    MAX_MOUNT_TRY=100

    if [ $# -eq 3 ]; then
      MOUNT_FILE=$3
      while [ ${MOUNT_READY} -eq 0 ]; do
        if mountpoint -q ${MOUNT_POINT} && [ -e ${MOUNT_POINT}/${MOUNT_FILE} ]; then
          MOUNT_READY=1
        else
          sleep 0.2
          let MOUNT_TRY++
          if [ $MOUNT_TRY -eq $MAX_MOUNT_TRY ]; then
	          MOUNT_READY=1
          fi  
        fi
      done
    else
      while [ ${MOUNT_READY} -eq 0 ]; do
        if mountpoint -q ${MOUNT_POINT}; then
          MOUNT_READY=1
        else
          sleep 0.2
          let MOUNT_TRY++
          if [ $MOUNT_TRY -eq $MAX_MOUNT_TRY ]; then
	          MOUNT_READY=1
          fi
        fi
      done  
    fi
  fi
}

if [ $(id -u) -ne 0 ]; then
  echo -e "\nPlease start the script as a root or sudo!\n"
  exit 1
else
  if [ $# -eq 2 ]; then

    FULL_PATH_TO_IMG_FILE=${1}
    MOUNT_POINT=${2}

    mkdir -p ${MOUNT_POINT}
    mount -o loop ${FULL_PATH_TO_IMG_FILE} ${MOUNT_POINT}
    sleep 1

    MOUNT_READY=0
    MOUNT_TRY=0
    MAX_MOUNT_TRY=100
    while [ ${MOUNT_READY} -eq 0 ]; do
      if mountpoint -q ${MOUNT_POINT} && [ -e ${MOUNT_POINT}/chroot-root/start.sh ]; then
          MOUNT_READY=1
      else
          sleep 0.2
          let MOUNT_TRY++
          if [ $MOUNT_TRY -eq $MAX_MOUNT_TRY ]; then
		      MOUNT_READY=1
			  exit 1
	      fi
      fi
    done

    if mountpoint -q ${MOUNT_POINT}; then
        mount_safe "${MOUNT_POINT}/chroot-root/gadgetron/proc" /proc self/exe
		mount_safe "${MOUNT_POINT}/chroot-root/gadgetron/dev" /dev
		mount_safe "${MOUNT_POINT}/chroot-root/gadgetron/sys" /sys
    else
      exit 1
    fi
    exit 0
  else
    echo -e "\nUsage: $0 <full path to img file> <mount point>\n"
    exit 1
  fi
fi
