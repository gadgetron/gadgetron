#!/bin/bash

image_name=$1
gadgetron_version=$(docker run --rm $image_name gadgetron --info | awk '/-- Version/ {print $4}')
manifest=$(docker run --rm $image_name /opt/code/manifest)
gadgetron_sha1=$(echo $manifest|jq '.io.gadgetron.gadgetron.sha1'| tr -d '"'|cut -c-8)
gtprep_sha1=$(echo $manifest|jq '.io.gadgetron.gtprep.sha1'| tr -d '"' | cut -c-8 )
tag_value="${gadgetron_version}-${gadgetron_sha1}"
if [ "$gtprep_sha1" != "null" ]; then
    tag_value="${tag_value}-${gtprep_sha1}"
fi
docker tag ${image_name} ${image_name}:${tag_value}
