#!/bin/bash

if [ $# -eq 2 ]; then

 if [ "$1" == "mount" ]; then
  sudo mount --bind /dev "${2}/dev"
  sudo mount --bind /sys "${2}/sys"
  sudo mount --bind /proc "${2}/proc"

  exit 0
 fi

 if [ "$1" == "umount" ]; then
  sudo umount "${2}/dev"
  sudo umount "${2}/sys"
  sudo umount "${2}/proc"

  exit 0
 fi

 echo -e "\nUsage:  $0 (mount or umount) (chrootdir)\n"
 exit 1

else
 echo -e "\nUsage:  $0 (mount or umount) (chrootdir)\n"
 exit 1
fi
