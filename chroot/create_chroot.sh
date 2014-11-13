#!/bin/bash

if [ $(id -u) -ne 0 ]; then
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1

else
 if [ $# -eq 5 ]; then

  # Add LIBRARY_PATHS to LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${5}

  export LC_ALL=C

  rm -rf ${2}/chroot/chroot-root

  mkdir -p ${2}/chroot/chroot-root

  touch ${2}/chroot/chroot-root/source-manifest.txt

  echo "gadgetron    ${4}" > ${2}/chroot/chroot-root/source-manifest.txt

  mkdir -p ${2}/chroot/chroot-root/gadgetron

  mkdir -p ${2}/chroot/chroot-root/gadgetron/webapp

  mkdir -p ${2}/chroot/chroot-backups

  apt-get install debootstrap -y

  debootstrap --variant=buildd --arch amd64 trusty ${2}/chroot/chroot-root/gadgetron http://gb.archive.ubuntu.com/ubuntu/

  cd ${2}
  make install DESTDIR="${2}/chroot/chroot-root/gadgetron" -j8

  ${3}/chroot/generate_gadgetron_root ${1} ${2}/chroot/chroot-root/gadgetron

  # Commented out for now. libcuda.so should be copied using copy-cuda-lib script.
  #cp ${6}/libcuda.so ${2}/chroot/chroot-root/gadgetron/${1}/lib  

  cp ${3}/chroot/start.sh ${2}/chroot/chroot-root
  cp ${3}/chroot/stop.sh ${2}/chroot/chroot-root
  cp ${3}/chroot/start-env.sh ${2}/chroot/chroot-root 
  cp ${3}/chroot/start-webapp.sh ${2}/chroot/chroot-root 
  cp ${3}/chroot/mount.sh ${2}/chroot/chroot-root

  chmod +x ${2}/chroot/copy-cuda-lib.sh
  cp ${2}/chroot/copy-cuda-lib.sh ${2}/chroot/chroot-root

  chmod +x ${2}/chroot/start-gadgetron.sh
  cp ${2}/chroot/start-gadgetron.sh ${2}/chroot/chroot-root/gadgetron 
 
  chmod +x ${2}/chroot/enter-chroot-env.sh
  cp ${2}/chroot/enter-chroot-env.sh ${2}/chroot/chroot-root/gadgetron

  chmod +x ${2}/chroot/run-webapp.sh
  cp ${2}/chroot/run-webapp.sh ${2}/chroot/chroot-root/gadgetron

  cp -n ${2}/chroot/chroot-root/gadgetron${1}/config/gadgetron.xml.example ${2}/chroot/chroot-root/gadgetron${1}/config/gadgetron.xml

  cp ${2}/chroot/gadgetron_web_app.cfg ${2}/chroot/chroot-root/gadgetron/webapp
  cp ${2}/chroot/gadgetron_web.conf ${2}/chroot/chroot-root/gadgetron/webapp
  cp ${2}/chroot/gadgetron_web_ld.conf ${2}/chroot/chroot-root/gadgetron/webapp
  cp ${3}/apps/gadgetron/webapp/gadgetron_web_app.py ${2}/chroot/chroot-root/gadgetron/webapp
  cp ${3}/apps/gadgetron/webapp/gadgetron_web_app.py ${2}/chroot/chroot-root/gadgetron${1}/bin/gadgetron_web_app.py
  cp ${2}/chroot/gadgetron_web_app.cfg ${2}/chroot/chroot-root/gadgetron${1}/config/

  chroot ${2}/chroot/chroot-root/gadgetron apt-get install python-dev python-twisted python-psutil -y 
  
  TAR_FILE_NAME=gadgetron-`date '+%Y%m%d-%H%M'`-${4:0:8}
  IMAGE_FILE_NAME=${2}/chroot/chroot-backups/${TAR_FILE_NAME}.img

  tar -zcf "${2}/chroot/chroot-backups/${TAR_FILE_NAME}.tar.gz" --directory "${2}/chroot" --exclude=./chroot-root/gadgetron/etc --exclude=./chroot-root/gadgetron/var --exclude=./chroot-root/gadgetron/dev --exclude=./chroot-root/gadgetron/sys --exclude=./chroot-root/gadgetron/proc --exclude=./chroot-root/gadgetron/root ./chroot-root

  dd if=/dev/zero of=${IMAGE_FILE_NAME} bs=1024k seek=1024 count=0
  mke2fs -F -t ext3 ${IMAGE_FILE_NAME}
  mkdir ${2}/chroot/gadgetron_root
  mount -o loop ${IMAGE_FILE_NAME} ${2}/chroot/gadgetron_root
  tar -xzf ${2}/chroot/chroot-backups/${TAR_FILE_NAME}.tar.gz -C ${2}/chroot/gadgetron_root/
  sleep 3
  umount ${2}/chroot/gadgetron_root
  rmdir ${2}/chroot/gadgetron_root

  rm -rf "${2}/chroot/chroot-root"

  chmod 666 "${2}/chroot/chroot-backups/${TAR_FILE_NAME}.tar.gz"
  chmod 666 "${IMAGE_FILE_NAME}"
 
  exit 0

 else
  echo -e "\nUsage:  $0 (gadgetron install prefix) (gadgetron binary dir) (gadgetron source dir) (GADGETRON_GIT_SHA1_HASH) (LIBRARY_PATHS) \n"
  exit 1
 fi

fi

# --ARGUMENTS--
# 1. cmake install prefix:     /usr/local/gadgetron
# 2. cmake binary dir:         /home/ubuntu/gadgetron/build
# 3. cmake source dir:         /home/ubuntu/gadgetron
# 4. HASH
# 5. LIBRARY_PATHS             /home/ubuntu/ismrmrd_install
