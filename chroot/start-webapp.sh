#!/bin/bash                                                                                                                                                    

#if [[ "$EUID" -ne 0 ]]; then
if [ $(id -u) -ne 0 ]; then 
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1

else
 BASEDIR=$(dirname $0)

 if [ $# -eq 1 ]; then

  mkdir -p "${1}/proc"
  if find "${1}/proc" -maxdepth 0 -empty | read v; then
   mount --bind /proc "${1}/proc";
  fi
  chroot "${1}" /run-webapp.sh
  exit 0

 elif [ $# -eq 0 ]; then

  mkdir -p $BASEDIR/gadgetron/proc
  if find $BASEDIR/gadgetron/proc -maxdepth 0 -empty | read v; then
   mount --bind /proc $BASEDIR/gadgetron/proc;
  fi
  chroot $BASEDIR/gadgetron /run-webapp.sh
  exit 0

 else
  echo -e "\nUsage: $0    or    $0 (chrootdir)\n"
  exit 1
 fi

fi
