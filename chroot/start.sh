#!/bin/bash

start_gadgetron_job=0
BASEDIR=$(dirname $0)

trap '(($start_gadgetron_job == 0)) || ((`kill -0 $start_gadgetron_job`)) || kill $start_gadgetron_job & while kill -0 $start_gadgetron_job 2>/dev/null; do sleep 1; done & $BASEDIR/stop.sh' HUP TERM INT

if [ $(id -u) -ne 0 ]; then
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1
else
 if [ $# -eq 0 ]; then
  $BASEDIR/mount.sh $BASEDIR
  chroot $BASEDIR/gadgetron /start-gadgetron.sh &
  start_gadgetron_job=($!)
  wait $!
  $BASEDIR/stop.sh
  exit 0
 else
  echo -e "\nUsage: $0\n"
  exit 1
 fi
fi
