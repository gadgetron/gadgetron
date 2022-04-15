#!/bin/bash

gadgetron_job=0
trap '(($gadgetron_job == 0)) || ((`kill -0 $gadgetron_job`))|| kill $gadgetron_job' HUP TERM INT

if [ $# -eq 1 ]; then
    LOG_FILE=${1}

    OMP_THREAD_LIMIT=$(nproc)
    PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/home/javeda2/anaconda3/envs/gadgetron/bin LD_LIBRARY_PATH=/home/javeda2/anaconda3/envs/gadgetron/lib:/usr/local/lib:/opt/intel/mkl/lib/intel64:/opt/intel/lib/intel64 /home/javeda2/anaconda3/envs/gadgetron/bin/gadgetron > ${LOG_FILE} &

    gadgetron_job=($!)
    wait $!
    exit 0
else
    echo -e "\nUsage: $0 (log file in chroot)\n"
    exit 1
fi
