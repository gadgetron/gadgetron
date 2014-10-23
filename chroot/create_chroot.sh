#!/bin/bash

if [ $# -eq 4 ]; then

 sudo rm -rf "${3}/chroot-root"

 mkdir -p "${3}/chroot-root"

 mkdir -p "${3}/chroot-root/gadgetron"

 mkdir -p "${3}/chroot-backups"

 sudo apt-get install debootstrap -y

 sudo debootstrap --variant=buildd --arch amd64 trusty "${3}/chroot-root/gadgetron" http://gb.archive.ubuntu.com/ubuntu/

 cd ${4}

 sudo make install DESTDIR="${3}/chroot-root/gadgetron" -j8

 "${2}/generate_gadgetron_root" ${1} "${3}/chroot-root/gadgetron"

 cp ${2}/start.sh "${3}/chroot-root"
 cp ${2}/stop.sh "${3}/chroot-root"
 cp ${2}/start-env.sh "${3}/chroot-root" 
 
 sudo chmod +x ${3}/start-gadgetron.sh
 sudo cp ${3}/start-gadgetron.sh "${3}/chroot-root/gadgetron" 
 
 sudo chmod +x ${3}/enter-chroot-env.sh
 sudo cp ${3}/enter-chroot-env.sh "${3}/chroot-root/gadgetron"

 sudo cp -n ${3}/chroot-root/gadgetron${1}/config/gadgetron.xml.example ${3}/chroot-root/gadgetron${1}/config/gadgetron.xml

 tar -zcf "${3}/chroot-backups/gadgetron-chroot-`date '+%d-%B-%Y'`.tar.gz" --directory "${3}/chroot-root" --exclude=./gadgetron/etc --exclude=./gadgetron/var --exclude=./gadgetron/dev --exclude=./gadgetron/root .

 sudo rm -rf "${3}/chroot-root"
 
 exit 0

else
 echo -e "\nUsage:  $0 (gadgetron install dir) (chroot source dir) (chroot build dir) (gadgetron build dir)\n"
 exit 1
fi
