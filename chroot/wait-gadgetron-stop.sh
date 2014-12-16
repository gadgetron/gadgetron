#!/bin/bash

function umount_check {
    MAX_TRY=100
    MOUNT_DIR=$1
    UMOUNT_READY=0
    UMOUNT_TRY=0
    while [ ${UMOUNT_READY} -eq 0 ]; do
        if mountpoint -q ${MOUNT_DIR}; then
            let UMOUNT_TRY++
            if [ $UMOUNT_TRY -eq $MAX_TRY ]; then
                UMOUNT_READY=1
            else
                sleep 0.2
            fi
        else
            UMOUNT_READY=1
        fi
    done
}

if [ $# -eq 1 ]; then

    MAX_TRY=100

    # check whether gadgetron has stopped
    CURRENT_TRY=0
    GADGETRON_PID=`pgrep -U root -x gadgetron`
    while [ -n "$GADGETRON_PID" ] && [ $GADGETRON_PID -ne 0 ]; do
        sleep 0.2
        GADGETRON_PID=`pgrep -U root -x gadgetron`

        if [ $CURRENT_TRY -eq $MAX_TRY ]; then
            GADGETRON_PID=0
        fi

        let CURRENT_TRY++
    done
     
    # check whether all umount actions has completed
    MOUNT_POINT=$1
    umount_check $MOUNT_POINT/chroot-root/gadgetron/proc
    umount_check $MOUNT_POINT/chroot-root/gadgetron/sys
    umount_check $MOUNT_POINT/chroot-root/gadgetron/dev
    umount_check $MOUNT_POINT

    exit 0

else
    echo -e "\nUsage: $0 (chrootdir)\n"
    exit 1
fi
