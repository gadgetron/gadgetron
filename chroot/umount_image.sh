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
}

if [ $(id -u) -ne 0 ]; then
  echo -e "\nPlease start the script as a root or sudo!\n"
  exit 1
else
  if [ $# -eq 1 ]; then

    MOUNT_POINT=${1}

    if mountpoint -q ${MOUNT_POINT}; then
      if mountpoint -q ${MOUNT_POINT}/chroot-root/gadgetron/proc; then
        umount ${MOUNT_POINT}/chroot-root/gadgetron/proc
        umount_check ${MOUNT_POINT}/chroot-root/gadgetron/proc
      fi
      if mountpoint -q ${MOUNT_POINT}/chroot-root/gadgetron/dev; then
        umount ${MOUNT_POINT}/chroot-root/gadgetron/dev
        umount_check ${MOUNT_POINT}/chroot-root/gadgetron/dev
      fi
      if mountpoint -q ${MOUNT_POINT}/chroot-root/gadgetron/sys; then
        umount ${MOUNT_POINT}/chroot-root/gadgetron/sys
        umount_check ${MOUNT_POINT}/chroot-root/gadgetron/sys
      fi
      if mountpoint -q ${MOUNT_POINT}/chroot-root/gadgetron/tmp/dependency; then
        umount ${MOUNT_POINT}/chroot-root/gadgetron/tmp/dependency
        umount_check ${MOUNT_POINT}/chroot-root/gadgetron/tmp/dependency
      fi
      umount ${MOUNT_POINT}
      umount_check ${MOUNT_POINT}
      exit 0
    fi
  else
    echo -e "\nUsage: $0 <mount point>\n"
    exit 1
  fi
fi
