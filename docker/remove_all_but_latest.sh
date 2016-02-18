#!/bin/bash

image_name=$1

docker rmi  $(docker images ${image_name} |grep -v "TAG"|grep -v "latest"|awk -v img=$image_name '{print img":"$2}')
