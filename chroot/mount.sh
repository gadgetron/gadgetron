#!/bin/bash                                                                                                                                    
if [ $(id -u) -ne 0 ]; then
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1

else
 BASEDIR=$(dirname $0)

 if [ $# -eq 1 ]; then

  mkdir -p "${1}/gadgetron/proc"
  if find "${1}/gadgetron/proc" -maxdepth 0 -empty | read v; then
   mount --bind /proc "${1}/gadgetron/proc";
  fi

  mkdir -p "${1}/gadgetron/dev"
  if find "${1}/gadgetron/dev" -maxdepth 0 -empty | read v; then
   mount --bind /dev "${1}/gadgetron/dev";
  fi

  mkdir -p "${1}/gadgetron/sys"
  if find "${1}/gadgetron/sys" -maxdepth 0 -empty | read v; then
   mount --bind /sys "${1}/gadgetron/sys";
  fi

  exit 0

 else
  echo -e "\nUsage: $0 (chrootdir)\n"
  exit 1
 fi

fi
