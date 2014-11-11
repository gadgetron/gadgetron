#!/bin/bash                                                                                                                                                  
if [ $(id -u) -ne 0 ]; then
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1

else
 BASEDIR=$(dirname $0)

 if [ $# -eq 0 ]; then
  $BASEDIR/mount.sh $BASEDIR
  chroot $BASEDIR/gadgetron /enter-chroot-env.sh
  exit 0

 else
  echo -e "\nUsage: $0\n"
  exit 1
 fi
fi
