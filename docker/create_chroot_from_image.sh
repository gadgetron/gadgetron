#!/bin/bash

image=$1
created=$(docker inspect --format='{{.Created}}' $image)
timestamp="${created:0:4}${created:5:2}${created:8:2}-${created:11:2}${created:14:2}"
gadgetron_hash=$(docker run --rm $image cut -c-8 /opt/code/gadgetron_sha1.txt)
gtprep_hash=$(docker run --rm $image bash -c  "[ -f /opt/code/gtprep_sha1.txt ] && cut -c-8 /opt/code/gtprep_sha1.txt")
chroot_name="gadgetron-${timestamp}-$gadgetron_hash"
if [ -n "$gtprep_hash" ]; then
    chroot_name="${chroot_name}-gtprep-${gtprep_hash}"
fi
chroot_base_name=$chroot_name
chroot_name="${chroot_name}.tar.gz"
echo "Creating chroot image with name: $chroot_name"

docker run -i --rm --volume=`pwd`:/opt/backup $image tar -cvpzf /opt/backup/$chroot_name --exclude=/opt/code --exclude=/opt/backup --one-file-system / 

if [ $# -gt 1 ]; then
    sleep 3
    chroot_img_name=${chroot_base_name}.img
    echo "Creating chroot img file with name: $chroot_img_name"

    chroot_img_size=$2
    chroot_tar_dir=${chroot_base_name}_untar
    mkdir -p ${chroot_tar_dir}
    tar -xzf ${chroot_name} --directory=${chroot_tar_dir}

    chroot_img_dir=${chroot_base_name}
    mkdir -p ${chroot_img_dir}
    dd if=/dev/zero of=${chroot_img_name} bs=${chroot_img_size}k seek=1024 count=0
    mke2fs -F -t ext3 ${chroot_img_name}
    mount -o loop ${chroot_img_name} ${chroot_img_dir}
    mkdir -p ${chroot_img_dir}/chroot-root
    mkdir -p ${chroot_img_dir}/chroot-root/gadgetron
    cp -rf ${chroot_tar_dir}/* ${chroot_img_dir}/chroot-root/gadgetron
    sleep 3
    umount ${chroot_img_dir}
    rmdir ${chroot_img_dir}
    rm -rf ${chroot_tar_dir}
fi
