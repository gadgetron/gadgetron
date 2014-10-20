#!/bin/bash

BASEDIR=$(dirname $0)

if [ $# -eq 1 ]; then

 sudo umount "${1}/proc"
 exit 0

elif [ $# -eq 0 ]; then

 sudo umount $BASEDIR/gadgetron/proc
 exit 0

else
 echo -e "\nUsage: $0    or    $0 (chrootdir)\n"
 exit 1
fi
