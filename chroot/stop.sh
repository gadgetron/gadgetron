#!/bin/bash

function umount_check {
    MAX_TRY=200
    MOUNT_DIR=$1
    UMOUNT_READY=0
    UMOUNT_TRY=0
    while [ ${UMOUNT_READY} -eq 0 ]; do
        if mountpoint -q ${MOUNT_DIR}; then
            let UMOUNT_TRY++
            if [ $UMOUNT_TRY -eq $MAX_TRY ]; then
                UMOUNT_READY=1
                umount $MOUNT_DIR
            else
                sleep 0.2
            fi
        else
            UMOUNT_READY=1
        fi
    done

    if mountpoint -q ${MOUNT_DIR}; then
        umount ${MOUNT_DIR}
    fi
}

if [ $(id -u) -ne 0 ]; then
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1

else
 BASEDIR=$(dirname $0)

 if [ $# -ge 1 ]; then
  CHROOT_DIR=${1} 

  if mountpoint -q $CHROOT_DIR/proc; then
   umount $CHROOT_DIR/proc
   umount_check $CHROOT_DIR/proc
  fi
  if mountpoint -q $CHROOT_DIR/sys; then
    umount $CHROOT_DIR/sys
   umount_check $CHROOT_DIR/sys
  fi
  if mountpoint -q $CHROOT_DIR/dev; then
    umount $CHROOT_DIR/dev
    umount_check $CHROOT_DIR/dev
  fi

  if [ $# -eq 2 ]; then
    DATA_DIR=${2}
    umount $CHROOT_DIR/tmp/gadgetron_data
    umount_check $CHROOT_DIR/tmp/gadgetron_data
  fi

  exit 0

 else
  echo -e "\nUsage: $0 (chroot_dir) (datadir)\n"
  exit 1
 fi
fi
