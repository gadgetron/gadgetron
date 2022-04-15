#!/bin/bash

if [ $(id -u) -ne 0 ]; then
 echo -e "\nPlease start the script as a root or sudo!\n"
 exit 1
else
 if [ $# -eq 1 ]; then
 BASEDIR=$(dirname $0)
 CHROOT_DIR=${1}

 # Absolute path this script is in
 SCRIPTPATH=$(dirname "$SCRIPT")

 # find the lib(s)
 CANDIDATES=$(ldconfig -p | grep "libcuda.so\s")

 # find the one that is 64-bit
 for CANDIDATE in $CANDIDATES 
 do
  var=$(file -L $CANDIDATE | grep '64-bit')
  if [ -n "$var" ]; then
   NEW_CUDA_LIB=$CANDIDATE
  fi 
 done

 # copy it to the right location (overwrite the previous one)
 yes | cp $NEW_CUDA_LIB $CHROOT_DIR/home/javeda2/anaconda3/envs/gadgetron/lib/
 exit 0

 else
  echo -e "\nUsage: $0 (chroot_dir)\n"
  exit 1
 fi
fi
