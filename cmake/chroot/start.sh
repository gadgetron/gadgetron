#!/bin/bash

start_gadgetron_job=0
BASEDIR=$(dirname $0)

trap '(($start_gadgetron_job == 0)) || ((`kill -0 $start_gadgetron_job`)) || kill $start_gadgetron_job & while kill -0 $start_gadgetron_job 2>/dev/null; do sleep 1; done & $BASEDIR/stop.sh $CHROOT_DIR' HUP TERM INT

if [ $(id -u) -ne 0 ]; then
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1
else
 if [ $# -eq 2 ]; then

  CHROOT_DIR=${1}
  LOG_FILE=${2}

  $BASEDIR/mount.sh $CHROOT_DIR
  chroot $CHROOT_DIR /home/javeda2/anaconda3/envs/gadgetron/share/gadgetron/chroot/start-gadgetron.sh ${LOG_FILE} &
  start_gadgetron_job=($!)
  wait $!
  $BASEDIR/stop.sh $CHROOT_DIR
  exit 0
 else
  echo -e "\nUsage: $0 (chroot_dir) (log file in chroot)\n"
  exit 1
 fi
fi
