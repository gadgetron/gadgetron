#!/bin/bash

if [ $(id -u) -ne 0 ]; then
  echo -e "\nPlease start the script as a root or sudo!\n"
  exit 1
else
  if [ $# -eq 0 ]; then
    if mountpoint -q /mnt/chroot; then
      if mountpoint -q /mnt/chroot/chroot-root/gadgetron/proc; then
        umount /mnt/chroot/chroot-root/gadgetron/proc
      fi
      if mountpoint -q /mnt/chroot/chroot-root/gadgetron/dev; then
        umount /mnt/chroot/chroot-root/gadgetron/dev
      fi
      if mountpoint -q /mnt/chroot/chroot-root/gadgetron/sys; then
        umount /mnt/chroot/chroot-root/gadgetron/sys
      fi
      umount /mnt/chroot
      exit 0
    fi
  else
    echo -e "\nUsage: $0\n"
    exit 1
  fi
fi
