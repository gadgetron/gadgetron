#!/bin/bash

if [ $(id -u) -ne 0 ]; then
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1

else
 if [ $# -eq 5 ]; then

  rm -rf ${3}/chroot-root

  mkdir -p ${3}/chroot-root

  mkdir -p ${3}/chroot-root/gadgetron

  mkdir -p ${3}/chroot-root/gadgetron/webapp

  mkdir -p ${3}/chroot-backups

  apt-get install debootstrap -y

  debootstrap --variant=buildd --arch amd64 trusty ${3}/chroot-root/gadgetron http://gb.archive.ubuntu.com/ubuntu/

  cd ${4}

  make install DESTDIR="${3}/chroot-root/gadgetron" -j8

  ${2}/generate_gadgetron_root ${1} ${3}/chroot-root/gadgetron

  cp ${2}/start.sh ${3}/chroot-root
  cp ${2}/stop.sh ${3}/chroot-root
  cp ${2}/start-env.sh ${3}/chroot-root 
  cp ${2}/start-webapp.sh ${3}/chroot-root 

  chmod +x ${3}/start-gadgetron.sh
  cp ${3}/start-gadgetron.sh ${3}/chroot-root/gadgetron 
 
  chmod +x ${3}/enter-chroot-env.sh
  cp ${3}/enter-chroot-env.sh ${3}/chroot-root/gadgetron

  chmod +x ${3}/run-webapp.sh
  cp ${3}/run-webapp.sh ${3}/chroot-root/gadgetron

  cp -n ${3}/chroot-root/gadgetron${1}/config/gadgetron.xml.example ${3}/chroot-root/gadgetron${1}/config/gadgetron.xml

  cp ${3}/gadgetron_web_app.cfg ${3}/chroot-root/gadgetron/webapp
  cp ${3}/gadgetron_web.conf ${3}/chroot-root/gadgetron/webapp
  cp ${3}/gadgetron_web_ld.conf ${3}/chroot-root/gadgetron/webapp
  cp ${5}/apps/gadgetron/webapp/gadgetron_web_app.py ${3}/chroot-root/gadgetron/webapp
  cp ${5}/apps/gadgetron/webapp/gadgetron_web_app.py ${3}/chroot-root/gadgetron${1}/bin/gadgetron_web_app.py
  cp ${3}/gadgetron_web_app.cfg ${3}/chroot-root/gadgetron${1}/config/

  chroot ${3}/chroot-root/gadgetron apt-get install python-dev python-twisted python-psutil -y 

  tar -zcf "${3}/chroot-backups/gadgetron-chroot-`date '+%d-%B-%Y'`.tar.gz" --directory "${3}" --exclude=./chroot-root/gadgetron/etc --exclude=./chroot-root/gadgetron/var --exclude=./chroot-root/gadgetron/dev --exclude=./chrot-root/gadgetron/root ./chroot-root

  TAR_FILE_NAME=gadgetron-chroot-${6}-`date '+%d-%B-%Y'`
  IMAGE_FILE_NAME=${3}/chroot-backups/${TAR_FILE_NAME}.img

  tar -zcf "${3}/chroot-backups/${TAR_FILE_NAME}.tar.gz" --directory "${3}" --exclude=./chroot-root/gadgetron/etc --exclude=./chroot-root/gadgetron/var --exclude=./chroot-root/gadgetron/dev --exclude=./chrot-root/gadgetron/root ./chroot-root

  dd if=/dev/zero of=${IMAGE_FILE_NAME} bs=1024k seek=1024 count=0
  mke2fs -F -t ext3 ${IMAGE_FILE_NAME}
  mkdir ${3}/gadgetron_root
  mount -o loop ${IMAGE_FILE_NAME} ${3}/gadgetron_root
  tar -xzf ${3}/chroot-backups/${TAR_FILE_NAME}.tar.gz -C ${3}/gadgetron_root/
  sleep 3
  umount ${3}/gadgetron_root
  rmdir ${3}/gadgetron_root

  rm -rf "${3}/chroot-root"
 
  exit 0

 else
  echo -e "\nUsage:  $0 (gadgetron install dir) (chroot source dir) (chroot build dir) (gadgetron build dir) (gadgetron source dir)\n"
  exit 1
 fi

fi
