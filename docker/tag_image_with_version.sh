#!/bin/bash

image_name=$1
gadgetron_version=$(docker run --rm $image_name gadgetron_info | awk '/-- Version/ {print $4}')
gadgetron_sha1=$(docker run --rm $image_name cut -c-8 /opt/code/gadgetron_sha1.txt)
docker tag ${image_name} ${image_name}:${gadgetron_version}-${gadgetron_sha1}
