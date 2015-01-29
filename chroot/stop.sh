#!/bin/bash

function umount_check {
    MAX_TRY=100
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
}

if [ $(id -u) -ne 0 ]; then
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1

else
 BASEDIR=$(dirname $0)

 if [ $# -eq 0 ]; then
  if mountpoint -q $BASEDIR/../proc; then
   umount $BASEDIR/../proc
   umount_check $BASEDIR/../proc
  fi
  if mountpoint -q $BASEDIR/../sys; then
    umount $BASEDIR/../sys
   umount_check $BASEDIR/../sys
  fi
  if mountpoint -q $BASEDIR/../dev; then
    umount $BASEDIR/../dev
    umount_check $BASEDIR/../dev
  fi
  exit 0

 else
  echo -e "\nUsage: $0\n"
  exit 1
 fi
fi
