#!/bin/bash

image_name=$1
gadgetron_version=$(docker run --rm $image_name gadgetron_info | awk '/-- Version/ {print $4}')
gadgetron_sha1=$(docker run --rm $image_name cut -c-8 /opt/code/gadgetron_sha1.txt)
gtprep_sha1=$(docker run --rm $image_name bash -c  "[ -f /opt/code/gtprep_sha1.txt ] && cut -c-8 /opt/code/gtprep_sha1.txt")
tag_value="${gadgetron_version}-${gadgetron_sha1}"
if [ -n "$gtprep_sha1" ]; then
    tag_value="${tag_value}-${gtprep_sha1}"
fi
docker tag ${image_name} ${image_name}:${tag_value}
