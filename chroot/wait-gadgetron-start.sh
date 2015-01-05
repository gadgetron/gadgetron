#!/bin/bash

BASEDIR=$(dirname $0)
  
CURRENT_TRY=0
MAX_TRY=100
GADGETRON_PID=`pgrep -U root -x gadgetron`
while [ -z "$GADGETRON_PID" ]; do
    sleep 0.2
    GADGETRON_PID=`pgrep -U root -x gadgetron`

    if [ $CURRENT_TRY -eq $MAX_TRY ]; then
        GADGETRON_PID=9999
    fi

    let CURRENT_TRY++
done

exit 0

