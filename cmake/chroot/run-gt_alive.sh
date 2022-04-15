#!/bin/bash

start_gadgetron_image_job=0
BASEDIR=$(dirname $0)

trap '(($start_gadgetron_image_job == 0)) || ((`kill -0 $start_gadgetron_image_job`)) || kill $start_gadgetron_image_job & while kill -0 $start_gadgetron_image_job 2>/dev/null; do sleep 1; done' HUP TERM INT

if [ $(id -u) -ne 0 ]; then
  echo -e "\nPlease start the script as a root or sudo!\n"
  exit 1
else
  if [ $# -ge 3 ]; then

    MOUNT_POINT=${1}
    HOSTNAME=${2}
    PORT=${3}

    if [ $# -eq 4 ]; then
        FULL_PATH_TO_IMG_FILE=${4}

        if find "${MOUNT_POINT}/chroot-root/gadgetron" -maxdepth 0 -empty | read v; then
            mkdir -p ${MOUNT_POINT}
            mount -o loop ${FULL_PATH_TO_IMG_FILE} ${MOUNT_POINT}
        fi
    fi

    chroot ${MOUNT_POINT}/chroot-root/gadgetron /home/javeda2/anaconda3/envs/gadgetron/share/gadgetron/chroot/gt_alive.sh $HOSTNAME $PORT &
	PID=$!
	(sleep 2 && kill -9 $PID) &
	waiter=$!
	wait $PID
	status=$?
	echo "gt_alive exit code : $status"
	kill -9 $waiter 2>/dev/null
	completeJob=$?
	
	if [ $completeJob -eq 0 ]; then
	    echo "gt_alive is completed properly"
	    exit $status
	else
	    exit 1
    fi
	
    exit $?
  else
    echo -e "\nUsage: $0 <mount point> <host> <port> <optional: full path to img file>\n"
    exit 1
  fi
fi
