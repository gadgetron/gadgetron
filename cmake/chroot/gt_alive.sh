#!/bin/bash

if [ $# -eq 2 ]; then

    GT_HOST=${1}
    GT_PORT=${2}

    PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/home/javeda2/anaconda3/envs/gadgetron/bin LD_LIBRARY_PATH=/home/javeda2/anaconda3/envs/gadgetron/lib:/usr/local/lib:/opt/intel/mkl/lib/intel64:/opt/intel/lib/intel64 /home/javeda2/anaconda3/envs/gadgetron/bin/gt_alive $GT_HOST $GT_PORT
    exit $?
else
    echo -e "\nUsage: $0 <host> <port>\n"
    exit 1
fi

exit 0
