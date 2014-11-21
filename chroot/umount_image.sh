#!/bin/bash

if [ $(id -u) -ne 0 ]; then
  echo -e "\nPlease start the script as a root or sudo!\n"
  exit 1
else
  if [ $# -eq 1 ]; then

    MOUNT_POINT=${1}

    if mountpoint -q ${MOUNT_POINT}; then
      if mountpoint -q ${MOUNT_POINT}/chroot-root/gadgetron/proc; then
        umount ${MOUNT_POINT}/chroot-root/gadgetron/proc
      fi
      if mountpoint -q ${MOUNT_POINT}/chroot-root/gadgetron/dev; then
        umount ${MOUNT_POINT}/chroot-root/gadgetron/dev
      fi
      if mountpoint -q ${MOUNT_POINT}/chroot-root/gadgetron/sys; then
        umount ${MOUNT_POINT}/chroot-root/gadgetron/sys
      fi
      umount ${MOUNT_POINT}
      exit 0
    fi
  else
    echo -e "\nUsage: $0 <mount point>\n"
    exit 1
  fi
fi
