#!/bin/bash

BASEDIR=$(dirname $0)

if [ $# -eq 0 ]; then
    CHROOT_ISMRMRD_DATA_PATH=/tmp/gadgetron_data
    CHROOT_ISMRMRD_DATA_HOURS=240
else
    if [ $# -eq 1 ]; then
        CHROOT_ISMRMRD_DATA_PATH=${1}
        CHROOT_ISMRMRD_DATA_HOURS=240
    else
        CHROOT_ISMRMRD_DATA_PATH=${1}
        CHROOT_ISMRMRD_DATA_HOURS=${2}
    fi
fi

# find all files in the ${CHROOT_ISMRMRD_DATA_PATH} folder

current_moment=`date +%s`

for file in ${CHROOT_ISMRMRD_DATA_PATH}/*.h5
do
    # check the file creation time, if this file is too old, delete it
    file_moment=`stat -c %Y ${file}`
    file_duration_seconds=$((current_moment-file_moment))
    file_duration_mins=$((file_duration_seconds / 60))
    file_duration_hours=$((file_duration_mins / 60))
    if [ ${file_duration_hours} -ge $CHROOT_ISMRMRD_DATA_HOURS ]; then
        rm -rf ${file}
    fi
done

exit 0
