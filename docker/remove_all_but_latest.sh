#!/bin/bash

image_name=$1

#Remove images not tagged with latest
rem_im=$(docker images ${image_name} |grep -v "TAG"|grep -v "latest"|awk -v img=$image_name '{print img":"$2}')
if [ -n "$rem_im" ]; then
    docker rmi  $rem_im
fi

#Remove any "orphaned" images
rem_im=$(docker images |grep "<none>"|awk  '{print $3}'|uniq)
if [ -n "$rem_im" ]; then
    docker rmi  $rem_im
fi

