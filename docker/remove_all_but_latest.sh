#!/bin/bash

image_name=$1

#Remove images not tagged with latest
docker rmi  $(docker images ${image_name} |grep -v "TAG"|grep -v "latest"|awk -v img=$image_name '{print img":"$2}')

#Remove any "orphaned" images
docker rmi  $(docker images |grep "<none>"|awk  '{print $3}'|uniq)
