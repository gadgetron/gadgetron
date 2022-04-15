#!/bin/bash

if [ $# -eq 5 ]; then

    ISMRMRD_FILENAME=${1}
    CONDIG_XML=${2}
    GT_HOST=${3}
    GT_PORT=${4}
    GT_TIMEOUT=${5}

    PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/home/javeda2/anaconda3/envs/gadgetron/bin LD_LIBRARY_PATH=/home/javeda2/anaconda3/envs/gadgetron/lib:/usr/local/lib:/opt/intel/mkl/lib/intel64:/opt/intel/lib/intel64 /home/javeda2/anaconda3/envs/gadgetron/bin/gadgetron_ismrmrd_client -f $ISMRMRD_FILENAME -c $CONDIG_XML -a $GT_HOST -p $GT_PORT -t $GT_TIMEOUT
    exit $?
else
    echo -e "\nUsage: $0 <ismrmrd filename> <config filename> <host> <port> <time out>\n"
    exit 1
fi

exit 0
