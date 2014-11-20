#!/bin/bash                                                                                                                                    
if [ $(id -u) -ne 0 ]; then
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1

else
 BASEDIR=$(dirname $0)

 if [ $# -eq 1 ]; then

  CHROOT_DIR=${1}

  mkdir -p "${CHROOT_DIR}/gadgetron/proc"
  if find "${CHROOT_DIR}/gadgetron/proc" -maxdepth 0 -empty | read v; then
   mount --bind /proc "${CHROOT_DIR}/gadgetron/proc";
  fi

  mkdir -p "${CHROOT_DIR}/gadgetron/dev"
  if find "${CHROOT_DIR}/gadgetron/dev" -maxdepth 0 -empty | read v; then
   mount --bind /dev "${CHROOT_DIR}/gadgetron/dev";
  fi

  mkdir -p "${CHROOT_DIR}/gadgetron/sys"
  if find "${CHROOT_DIR}/gadgetron/sys" -maxdepth 0 -empty | read v; then
   mount --bind /sys "${CHROOT_DIR}/gadgetron/sys";
  fi

  exit 0

 else
  echo -e "\nUsage: $0 (chrootdir)\n"
  exit 1
 fi

fi
