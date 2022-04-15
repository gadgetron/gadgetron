#!/bin/bash

if [ $# -eq 3 ]; then

    DAT_FILENAME=${1}
    ISMRMRD_FILENAME=${2}
    SCAN_NO=${3}
    DATE=`date +%Y-%m-%d`

    PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/home/javeda2/anaconda3/envs/gadgetron/bin LD_LIBRARY_PATH=/home/javeda2/anaconda3/envs/gadgetron/lib:/usr/local/lib:/opt/intel/mkl/lib/intel64:/opt/intel/lib/intel64 /usr/local/bin/siemens_to_ismrmrd -f $DAT_FILENAME -o $ISMRMRD_FILENAME -z $SCAN_NO --studyDate $DATE
    exit $?
else
    echo -e "\nUsage: $0 <dat filename> <ismrmrd filename> <scan number>\n"
    exit 1
fi

exit 0
