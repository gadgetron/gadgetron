#!/bin/bash

BASEDIR=$(dirname $0)

if [ $# -eq 0 ]; then
    CHROOT_ISMRMRD_DATA_PATH=/tmp/gadgetron_data
    CHROOT_ISMRMRD_DISK_USAGE=80
    CHROOT_ISMRMRD_DATA_HOURS=72
else
    if [ $# -eq 1 ]; then
        CHROOT_ISMRMRD_DATA_PATH=${1}
        CHROOT_ISMRMRD_DISK_USAGE=80
        CHROOT_ISMRMRD_DATA_HOURS=72
    else
        if [ $# -eq 2 ]; then
            CHROOT_ISMRMRD_DATA_PATH=${1}
            CHROOT_ISMRMRD_DISK_USAGE=${2}
            CHROOT_ISMRMRD_DATA_HOURS=72
        else
            CHROOT_ISMRMRD_DATA_PATH=${1}
            CHROOT_ISMRMRD_DISK_USAGE=${2}
            CHROOT_ISMRMRD_DATA_HOURS=${3}
        fi
    fi
fi

# check disk usage
disk_used=`df /tmp/gadgetron_data | awk '$1!="Filesystem"{print $3}'`
disk_available=`df /tmp/gadgetron_data | awk '$1!="Filesystem"{print $4}'`
disk_total=$((disk_available+disk_used))
disk_usage_ratio=$((100*disk_used/disk_total + 1))
echo Disk usage ratio: ${disk_usage_ratio}%
echo --------------------------------------

if [ $disk_usage_ratio -ge $CHROOT_ISMRMRD_DISK_USAGE ]; then

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
            echo Delete --- ${file}
        else
            echo Keep  --- ${file}
        fi
    done
else
    for file in ${CHROOT_ISMRMRD_DATA_PATH}/*.h5
    do
        echo Keep  --- ${file}
    done
fi
echo --------------------------------------

exit 0
