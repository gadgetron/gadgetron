#!/bin/bash

start_gadgetron_image_job=0
BASEDIR=$(dirname $0)

trap '(($start_gadgetron_image_job == 0)) || ((`kill -0 $start_gadgetron_image_job`)) || kill $start_gadgetron_image_job & while kill -0 $start_gadgetron_image_job 2>/dev/null; do sleep 1; done & $BASEDIR/umount_image.sh' HUP TERM INT

if [ $(id -u) -ne 0 ]; then
  echo -e "\nPlease start the script as a root or sudo!\n"
  exit 1
else
  if [ $# -eq 2 ]; then
    mkdir -p ${2}
    mount -o loop ${1} ${2}
    ${2}/chroot-root/start.sh &
    start_gadgetron_image_job=($!)
    wait $!
    sleep 1
    $BASEDIR/umount_image.sh
    exit 0
  else
    echo -e "\nUsage: $0 <full path to img file> <mount point>\n"
    exit 1
  fi
fi
