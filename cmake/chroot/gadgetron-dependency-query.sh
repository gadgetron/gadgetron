#!/bin/bash

if [ $# -eq 4 ]; then

    GT_HOST=${1}
    GT_PORT=${2}
    QUERY_OUT=${3}
    TIMEOUT=${4}

    PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/home/javeda2/anaconda3/envs/gadgetron/bin LD_LIBRARY_PATH=/home/javeda2/anaconda3/envs/gadgetron/lib:/usr/local/lib:/opt/intel/mkl/lib/intel64:/opt/intel/lib/intel64 /home/javeda2/anaconda3/envs/gadgetron/bin/gadgetron_ismrmrd_client -q -c gtquery.xml -a $GT_HOST -p $GT_PORT -t $TIMEOUT -o $QUERY_OUT
    exit $?
else
    echo -e "\nUsage: $0 <Host> <port> <query out file> <time out in ms>\n"
    exit 1
fi

exit 0
