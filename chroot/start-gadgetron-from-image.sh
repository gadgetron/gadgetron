#!/bin/bash

start_gadgetron_image_job=0
BASEDIR=$(dirname $0)

trap '(($start_gadgetron_image_job == 0)) || ((`kill -0 $start_gadgetron_image_job`)) || kill $start_gadgetron_image_job & while kill -0 $start_gadgetron_image_job 2>/dev/null; do sleep 1; done' HUP TERM INT

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

    if mountpoint -q ${MOUNT_POINT}; then
	  ${MOUNT_POINT}/chroot-root/start.sh &
	  start_gadgetron_image_job=($!)
	  wait $!
	  $BASEDIR/umount_image.sh ${MOUNT_POINT}
    else
      exit 1
    fi
    exit 0
  else
    echo -e "\nUsage: $0 <full path to img file> <mount point>\n"
    exit 1
  fi
fi
