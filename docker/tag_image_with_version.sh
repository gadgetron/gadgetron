#!/bin/bash

image_name=$1
gadgetron_version=$(docker run --rm $image_name gadgetron_info | awk '/-- Version/ {print $4}')
docker tag ${image_name} ${image_name}:${gadgetron_version}
