#!/bin/bash

if [ $(id -u) -ne 0 ]; then
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1

else
 if [ $# -eq 6 ]; then

# --ARGUMENTS--
# 1. CHROOT_GADGETRON_INSTALL_PREFIX:     /usr/local/gadgetron
# 2. CHROOT_GADGETRON_BINARY_DIR:         /home/ubuntu/gadgetron/build
# 3. CHROOT_GADGETRON_SOURCE_DIR:         /home/ubuntu/gadgetron
# 4. CHROOT_GIT_SHA1_HASH
# 5. CHROOT_LIBRARY_PATHS
# 6. CHROOT_CUDA_LIBRARY

  CHROOT_GADGETRON_INSTALL_PREFIX=${1}
  echo CHROOT_GADGETRON_INSTALL_PREFIX: ${CHROOT_GADGETRON_INSTALL_PREFIX}
  
  CHROOT_GADGETRON_BINARY_DIR=${2}
  echo CHROOT_GADGETRON_BINARY_DIR: ${CHROOT_GADGETRON_BINARY_DIR}
  
  CHROOT_GADGETRON_SOURCE_DIR=${3}
  echo CHROOT_GADGETRON_SOURCE_DIR: ${CHROOT_GADGETRON_SOURCE_DIR}
  
  CHROOT_GIT_SHA1_HASH=${4}
  echo CHROOT_GIT_SHA1_HASH: ${CHROOT_GIT_SHA1_HASH}
  
  CHROOT_LIBRARY_PATHS=${5}
  echo CHROOT_LIBRARY_PATHS: ${CHROOT_LIBRARY_PATHS}
  
  CHROOT_CUDA_LIBRARY=${6}
  echo CHROOT_CUDA_LIBRARY: ${CHROOT_CUDA_LIBRARY}

  # Add LIBRARY_PATHS to LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CHROOT_LIBRARY_PATHS}
  export LC_ALL=C

  rm -rf ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root

  mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root

  touch ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/source-manifest.txt

  echo "gadgetron    ${CHROOT_GIT_SHA1_HASH}" > ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/source-manifest.txt

  mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron

  mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/webapp

  mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups

  apt-get install debootstrap -y

  debootstrap --variant=buildd --arch amd64 trusty ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron http://gb.archive.ubuntu.com/ubuntu/

  cd ${CHROOT_GADGETRON_BINARY_DIR}
  make install DESTDIR="${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron" -j8

  ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/generate_gadgetron_root ${CHROOT_GADGETRON_INSTALL_PREFIX} ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron

  cp ${6} ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/${CHROOT_GADGETRON_INSTALL_PREFIX}/lib  

  cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/start.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root
  cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/stop.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root
  cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/start-env.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root 
  cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/start-webapp.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root 
  cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/mount.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root

  cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/umount_image.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups
  cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/start-gadgetron-from-image.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups

  chmod +x ${CHROOT_GADGETRON_BINARY_DIR}/chroot/copy-cuda-lib.sh
  cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/copy-cuda-lib.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root

  chmod +x ${CHROOT_GADGETRON_BINARY_DIR}/chroot/start-gadgetron.sh
  cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/start-gadgetron.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron 
 
  chmod +x ${CHROOT_GADGETRON_BINARY_DIR}/chroot/enter-chroot-env.sh
  cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/enter-chroot-env.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron

  chmod +x ${CHROOT_GADGETRON_BINARY_DIR}/chroot/run-webapp.sh
  cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/run-webapp.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron

  cp -n ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron${CHROOT_GADGETRON_INSTALL_PREFIX}/config/gadgetron.xml.example ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron${CHROOT_GADGETRON_INSTALL_PREFIX}/config/gadgetron.xml

  cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_web_app.cfg ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/webapp
  cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_web.conf ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/webapp
  cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_web_ld.conf ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/webapp
  cp ${CHROOT_GADGETRON_SOURCE_DIR}/apps/gadgetron/webapp/gadgetron_web_app.py ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/webapp
  cp ${CHROOT_GADGETRON_SOURCE_DIR}/apps/gadgetron/webapp/gadgetron_web_app.py ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron${CHROOT_GADGETRON_INSTALL_PREFIX}/bin/gadgetron_web_app.py
  cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_web_app.cfg ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron${CHROOT_GADGETRON_INSTALL_PREFIX}/config/

  chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get install python-dev python-twisted python-psutil -y 
  
  TAR_FILE_NAME=gadgetron-`date '+%Y%m%d-%H%M'`-${CHROOT_GIT_SHA1_HASH:0:8}
  IMAGE_FILE_NAME=${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/${TAR_FILE_NAME}.img

  tar -zcf "${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/${TAR_FILE_NAME}.tar.gz" --directory "${CHROOT_GADGETRON_BINARY_DIR}/chroot" --exclude=./chroot-root/gadgetron/etc --exclude=./chroot-root/gadgetron/var --exclude=./chroot-root/gadgetron/dev --exclude=./chroot-root/gadgetron/sys --exclude=./chroot-root/gadgetron/proc --exclude=./chroot-root/gadgetron/root ./chroot-root

  dd if=/dev/zero of=${IMAGE_FILE_NAME} bs=1024k seek=1024 count=0
  mke2fs -F -t ext3 ${IMAGE_FILE_NAME}
  mkdir ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_root
  mount -o loop ${IMAGE_FILE_NAME} ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_root
  tar -xzf ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/${TAR_FILE_NAME}.tar.gz -C ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_root/
  sleep 3
  umount ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_root
  rmdir ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_root

  rm -rf "${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root"

  chmod 666 "${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/${TAR_FILE_NAME}.tar.gz"
  chmod 666 "${IMAGE_FILE_NAME}"
 
  exit 0

 else
  echo -e "\nUsage:  $0 (gadgetron install prefix) (gadgetron binary dir) (gadgetron source dir) (GADGETRON_GIT_SHA1_HASH) (LIBRARY_PATHS) (CUDA_LIBRARY)\n"
  exit 1
 fi

fi
