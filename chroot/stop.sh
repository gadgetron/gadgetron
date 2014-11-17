#!/bin/bash                                                                                                                                                 

if [ $(id -u) -ne 0 ]; then
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1

else
 BASEDIR=$(dirname $0)

 if [ $# -eq 0 ]; then
  if mountpoint -q $BASEDIR/gadgetron/proc; then
   umount $BASEDIR/gadgetron/proc
  fi
  if mountpoint -q $BASEDIR/gadgetron/sys; then
    umount $BASEDIR/gadgetron/sys
  fi
  if mountpoint -q $BASEDIR/gadgetron/dev; then
    umount $BASEDIR/gadgetron/dev
  fi
  exit 0

 else
  echo -e "\nUsage: $0\n"
  exit 1
 fi
fi
