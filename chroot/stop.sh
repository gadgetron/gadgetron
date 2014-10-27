#!/bin/bash

if [[ "$EUID" -ne 0 ]]; then
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1

else
 BASEDIR=$(dirname $0)

 if [ $# -eq 1 ]; then

  umount "${1}/proc"
  exit 0

 elif [ $# -eq 0 ]; then

  umount $BASEDIR/gadgetron/proc
  exit 0

 else
  echo -e "\nUsage: $0    or    $0 (chrootdir)\n"
  exit 1
 fi

fi
