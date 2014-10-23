#!/bin/bash                                                                                                                                                  

BASEDIR=$(dirname $0)

if [ $# -eq 1 ]; then

 sudo mkdir -p "${1}/proc"
 if find "${1}/proc" -maxdepth 0 -empty | read v; then
  sudo mount --bind /proc "${1}/proc";
 fi
 sudo chroot "${1}" /enter-chroot-env.sh
 exit 0

elif [ $# -eq 0 ]; then

 sudo mkdir -p $BASEDIR/gadgetron/proc
 if find $BASEDIR/gadgetron/proc -maxdepth 0 -empty | read v; then
  sudo mount --bind /proc $BASEDIR/gadgetron/proc;
 fi
 sudo chroot $BASEDIR/gadgetron /enter-chroot-env.sh
 exit 0

else
 echo -e "\nUsage: $0    or    $0 (chrootdir)\n"
 exit 1
fi
