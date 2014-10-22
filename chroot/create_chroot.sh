#!/bin/bash

if [ $# -eq 3 ]; then

 cd ${3} 

 sudo rm -rf "${2}/chroot-root"

 mkdir -p "${2}/chroot-root"

 mkdir -p "${2}/chroot-root/gadgetron"

 mkdir -p "${2}/chroot-backups"

 sudo apt-get install debootstrap -y

 sudo debootstrap --variant=buildd --arch amd64 trusty "${2}/chroot-root/gadgetron" http://gb.archive.ubuntu.com/ubuntu/

 sudo make install DESTDIR="${2}/chroot-root/gadgetron" -j8

 "${2}/generate_gadgetron_root" ${1} "${2}/chroot-root/gadgetron"

 cp ${2}/start.sh "${2}/chroot-root"
 cp ${2}/stop.sh "${2}/chroot-root"
 sudo cp ${2}/start-gadgetron.sh "${2}/chroot-root/gadgetron" 

 sudo cp -n ${2}/chroot-root/gadgetron${1}/config/gadgetron.xml.example ${2}/chroot-root/gadgetron${1}/config/gadgetron.xml

 tar -zcf "${2}/chroot-backups/gadgetron-chroot-`date '+%d-%B-%Y'`.tar.gz" --directory "${2}/chroot-root" --exclude=./gadgetron/etc --exclude=./gadgetron/var --exclude=./gadgetron/dev --exclude=./gadgetron/root .

 sudo rm -rf "${2}/chroot-root"
 
 exit 0

else
 echo -e "\nUsage:  $0 (gadgetron_root) (chroot scripts dir) (gadgetron build dir)\n"
 exit 1
fi
