#!/bin/bash

if [ $(id -u) -ne 0 ]; then 
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1

else
 BASEDIR=$(dirname $0)

 if [ $# -eq 0 ]; then
  umount $BASEDIR/gadgetron/proc
  umount $BASEDIR/gadgetron/sys
  umount $BASEDIR/gadgetron/dev
  exit 0

 else
  echo -e "\nUsage: $0\n"
  exit 1
 fi
fi
