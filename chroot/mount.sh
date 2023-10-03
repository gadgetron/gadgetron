#!/bin/bash                                                                                                                                    

function mount_safe {
  MOUNT_POINT=$1
  MOUNT_DIR=$2
  mkdir -p $MOUNT_POINT
  if find $MOUNT_POINT -maxdepth 0 -empty | read v; then
    mount --bind $MOUNT_DIR $MOUNT_POINT
    MOUNT_READY=0
    MOUNT_TRY=0
    MAX_MOUNT_TRY=100

    if [ $# -eq 3 ]; then
      MOUNT_FILE=$3
      while [ ${MOUNT_READY} -eq 0 ]; do
        if mountpoint -q ${MOUNT_POINT} && [ -e ${MOUNT_POINT}/${MOUNT_FILE} ]; then
          MOUNT_READY=1
        else
          sleep 0.2
          let MOUNT_TRY++
          if [ $MOUNT_TRY -eq $MAX_MOUNT_TRY ]; then
            MOUNT_READY=1
          fi  
        fi
      done
    else
      while [ ${MOUNT_READY} -eq 0 ]; do
        if mountpoint -q ${MOUNT_POINT}; then
          MOUNT_READY=1
        else
          sleep 0.2
          let MOUNT_TRY++
          if [ $MOUNT_TRY -eq $MAX_MOUNT_TRY ]; then
            MOUNT_READY=1
          fi
        fi
      done  
    fi
  fi
}

if [ $(id -u) -ne 0 ]; then
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1

else
 BASEDIR=$(dirname $0)

 if [ $# -ge 1 ]; then

  CHROOT_DIR=${1}

  mount_safe "${CHROOT_DIR}/proc" /proc self/exe
  mount_safe "${CHROOT_DIR}/dev" /dev
  mount_safe "${CHROOT_DIR}/sys" /sys

  if [ $# -eq 2 ]; then
    DATA_DIR=${2}
    mkdir -p ${CHROOT_DIR}/tmp/gadgetron_data
    mount_safe "${CHROOT_DIR}/tmp/gadgetron_data" ${DATA_DIR}
  fi

  exit 0

 else
  echo -e "\nUsage: $0 (chrootdir) (datadir)\n"
  exit 1
 fi

fi
