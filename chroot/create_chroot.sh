#!/bin/bash

if [ $(id -u) -ne 0 ]; then
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1

else
 if [ $# -ge 6 ]; then

# --ARGUMENTS-- (example)

# CHROOT_GADGETRON_INSTALL_PREFIX:    /usr/local/gadgetron
# CHROOT_GADGETRON_BINARY_DIR:        /home/ubuntu/gadgetron/build
# CHROOT_GADGETRON_SOURCE_DIR:        /home/ubuntu/gadgetron
# CHROOT_GIT_SHA1_HASH:               f4d7a9189fd21b07e482d28ecb8b07e589f81f9e
# CHROOT_LIBRARY_PATHS:               /usr/local/lib:/usr/lib/x86_64-linux-gnu
# CHROOT_CUDA_LIBRARY:                /usr/lib/x86_64-linux-gnu/libcuda.so

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

  echo "***** LD_LIBRARY_PATH ***** : ${LD_LIBRARY_PATH}"

  rm -rf ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root

  mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root

  touch ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/source-manifest.txt

  echo "gadgetron    ${CHROOT_GIT_SHA1_HASH}" > ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/source-manifest.txt

  mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron

  mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/webapp

  mkdir -p ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups

  apt-get install debootstrap -y

  debootstrap --variant=buildd --arch amd64 trusty ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron http://gb.archive.ubuntu.com/ubuntu/

  chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get --yes install software-properties-common python-dev python-twisted python-psutil python-numpy

  chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron add-apt-repository "deb http://us-east-1.ec2.archive.ubuntu.com/ubuntu/ trusty restricted main multiverse universe"
  chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron add-apt-repository "deb http://us-east-1.ec2.archive.ubuntu.com/ubuntu/ trusty-updates universe restricted multiverse main"
  chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron add-apt-repository "deb http://security.ubuntu.com/ubuntu trusty-security main universe"
  chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron add-apt-repository "deb http://gb.archive.ubuntu.com/ubuntu trusty main"
  chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron add-apt-repository "deb http://gadgetronubuntu.s3.amazonaws.com trusty main"

  chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get update
  chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get install libopenblas-base
  chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get --yes install sudo
  chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get --yes install python-h5py

  chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get --yes --force-yes install siemens-to-ismrmrd
  chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get --yes --force-yes install gadgetron-whole
  chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get --yes --force-yes install gadgetron-scripts

  #cp ${6} ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/${CHROOT_GADGETRON_INSTALL_PREFIX}/lib  
  #cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/start.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root
  #cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/stop.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root
  #cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/start-env.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root 
  #cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/start-webapp.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root 
  #cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/mount.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root
  #cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/mount_image.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups
  #cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/umount_image.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups
  #cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/start-gadgetron-from-image.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups
  #cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/run-gadgetron_ismrmrd_client.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups
  #cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/run-siemens_to_ismrmrd.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups
  #cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/run-gadgetron-dependency-query.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups
  #cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/run-gt_alive.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups

  #chmod +x ${CHROOT_GADGETRON_BINARY_DIR}/chroot/copy-cuda-lib.sh
  #cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/copy-cuda-lib.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root
  #chmod +x ${CHROOT_GADGETRON_BINARY_DIR}/chroot/start-gadgetron.sh
  #cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/start-gadgetron.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron 
  #chmod +x ${CHROOT_GADGETRON_BINARY_DIR}/chroot/enter-chroot-env.sh
  #cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/enter-chroot-env.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron
  #chmod +x ${CHROOT_GADGETRON_BINARY_DIR}/chroot/run-webapp.sh
  #cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/run-webapp.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron
  #chmod +x ${CHROOT_GADGETRON_BINARY_DIR}/chroot/siemens_to_ismrmrd.sh
  #cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/siemens_to_ismrmrd.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron
  #chmod +x ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gt_alive.sh
  #cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gt_alive.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron
  #chmod +x ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_ismrmrd_client.sh
  #cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_ismrmrd_client.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron
  #chmod +x ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron-dependency-query.sh
  #cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron-dependency-query.sh ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron
  cp -n ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron${CHROOT_GADGETRON_INSTALL_PREFIX}/config/gadgetron.xml.example ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron${CHROOT_GADGETRON_INSTALL_PREFIX}/config/gadgetron.xml
  cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_web_app.cfg ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/webapp
  cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_web.conf ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/webapp
  cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_web_ld.conf ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/webapp
  cp ${CHROOT_GADGETRON_SOURCE_DIR}/chroot/gadgetron_chroot.conf ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/webapp/gadgetron_chroot.conf
  cp ${CHROOT_GADGETRON_SOURCE_DIR}/apps/gadgetron/webapp/gadgetron_web_app.py ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron/webapp
  cp ${CHROOT_GADGETRON_SOURCE_DIR}/apps/gadgetron/webapp/gadgetron_web_app.py ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron${CHROOT_GADGETRON_INSTALL_PREFIX}/bin/gadgetron_web_app.py
  cp ${CHROOT_GADGETRON_BINARY_DIR}/chroot/gadgetron_web_app.cfg ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron${CHROOT_GADGETRON_INSTALL_PREFIX}/config/

  chroot ${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-root/gadgetron apt-get clean

  TAR_FILE_NAME=gadgetron-`date '+%Y%m%d-%H%M'`-${CHROOT_GIT_SHA1_HASH:0:8}
  IMAGE_FILE_NAME=${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/${TAR_FILE_NAME}.img

  tar -zcf "${CHROOT_GADGETRON_BINARY_DIR}/chroot/chroot-backups/${TAR_FILE_NAME}.tar.gz" --directory "${CHROOT_GADGETRON_BINARY_DIR}/chroot" --exclude=./chroot-root/gadgetron/dev --exclude=./chroot-root/gadgetron/sys --exclude=./chroot-root/gadgetron/proc --exclude=./chroot-root/gadgetron/root ./chroot-root

  dd if=/dev/zero of=${IMAGE_FILE_NAME} bs=1536k seek=1024 count=0
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
