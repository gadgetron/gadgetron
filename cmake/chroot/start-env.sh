#!/bin/bash                                                                                                                                                  
if [ $(id -u) -ne 0 ]; then
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1

else
 BASEDIR=$(dirname $0)

 if [ $# -eq 1 ]; then

  CHROOT_DIR=${1}
  
  $BASEDIR/mount.sh $CHROOT_DIR
  chroot $CHROOT_DIR/ /home/javeda2/anaconda3/envs/gadgetron/share/gadgetron/chroot/enter-chroot-env.sh
  exit 0

 else
  echo -e "\nUsage: $0 (chroot dir)\n"
  exit 1
 fi
fi
