#!/bin/bash

if [ $(id -u) -ne 0 ]; then
  echo -e "\nPlease start the script as a root or sudo!\n"
  exit 1
else
  if [ $# -eq 1 ]; then
    if mountpoint -q ${1}; then
      if mountpoint -q ${1}/chroot-root/gadgetron/proc; then
        umount ${1}/chroot-root/gadgetron/proc
      fi
      if mountpoint -q ${1}/chroot-root/gadgetron/dev; then
        umount ${1}/chroot-root/gadgetron/dev
      fi
      if mountpoint -q ${1}/chroot-root/gadgetron/sys; then
        umount ${1}/chroot-root/gadgetron/sys
      fi
      umount ${1}
      exit 0
    fi
  else
    echo -e "\nUsage: $0 <mount point>\n"
    exit 1
  fi
fi
