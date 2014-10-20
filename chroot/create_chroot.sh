#!/bin/bash

if [ $# -eq 3 ]; then

 rm -rf "${2}"

 mkdir -p "${2}"

 mkdir -p "${3}/chroot-backups"

 apt-get install debootstrap -y

 debootstrap --variant=buildd --arch amd64 trusty "${2}" http://gb.archive.ubuntu.com/ubuntu/

 "${3}/generate_gadgetron_root" "${1}" "${2}"

 tar -zcvpf "${3}/chroot-backups/full-backup-`date '+%d-%B-%Y'`.tar.gz" --directory "${2}" --exclude="${2}/proc" --exclude="${2}/sys" --exclude="${2}/dev" .

 rm -rf "${2}"
 
 exit 0

else
 echo -e "\nUsage:  $0 (gadgetron_root) (chroot dir) (chroot scripts dir)\n"
 exit 1
fi
