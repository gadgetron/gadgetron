#!/bin/bash

image=$1
created=$(docker inspect --format='{{.Created}}' $image)
timestamp="${created:0:4}${created:5:2}${created:8:2}-${created:11:2}${created:14:2}"
gadgetron_hash=$(docker run --rm $image cut -c-8 /opt/code/gadgetron_sha1.txt)
gtprep_hash=$(docker run --rm $image [ -f /opt/code/gtprep_sha1.txt ] && $(cut -c-8 /opt/code/gtprep_sha1.txt))
chroot_name="gadgetron-${timestamp}-$gadgetron_hash"
if [ -n "$gtprep_hash" ]; then
    chroot_name="${chroot_name}-gtprep-${gtprep_hash}"
fi
chroot_name="${chroot_name}.tar.gz"
echo "Creating chroot image with name: $chroot_name"

docker run -ti --rm --volume=`pwd`:/opt/backup $image tar -cvpzf /opt/backup/$chroot_name --exclude=/opt/code --exclude=/opt/backup --one-file-system / 
