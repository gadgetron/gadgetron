#!/bin/sh

NVIDIA_DEVICES=""

for f in `ls -a /dev/* | grep nvidia`
do
    NVIDIA_DEVICES="${NVIDIA_DEVICES} --device $f:$f"
done

echo $NVIDIA_DEVICES
