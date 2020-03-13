#!/bin/bash

base_dir=$(pwd)

cd ${base_dir}/base/ubuntu_1604
docker build --network=host -t gadgetron/ubuntu_1604_base .

cd ${base_dir}/base/ubuntu_1604_basic
docker build --network=host -t gadgetron/ubuntu_1604_base_basic .

#cd ${base_dir}/base/ubuntu_1604_cuda90_cudnn7
#docker build --network=host -t gadgetron/ubuntu_1604_cuda90_cudnn7_base .

#cd ${base_dir}/base/ubuntu_1604_cuda92_cudnn7
#docker build --network=host -t gadgetron/ubuntu_1604_cuda92_cudnn7_base .

cd ${base_dir}/base/ubuntu_1604_cuda100_cudnn7
docker build --network=host -t gadgetron/ubuntu_1604_cuda100_cudnn7_base .

cd ${base_dir}/base/ubuntu_1604_cuda101_cudnn7
docker build --network=host -t gadgetron/ubuntu_1604_cuda101_cudnn7_base .

cd ${base_dir}/base/ubuntu_1804
docker build --network=host -t gadgetron/ubuntu_1804_base .

cd ${base_dir}/base/ubuntu_1804_basic
docker build --network=host -t gadgetron/ubuntu_1804_base_basic .

#cd ${base_dir}/base/ubuntu_1804_cuda90
#docker build --network=host -t gadgetron/ubuntu_1804_cuda90_cudnn7_base .

#cd ${base_dir}/base/ubuntu_1804_cuda92_cudnn7
#docker build --network=host -t gadgetron/ubuntu_1804_cuda92_cudnn7_base .

cd ${base_dir}/base/ubuntu_1804_cuda100_cudnn7
docker build --network=host -t gadgetron/ubuntu_1804_cuda100_cudnn7_base .

cd ${base_dir}/base/ubuntu_1804_cuda101_cudnn7
docker build --network=host -t gadgetron/ubuntu_1804_cuda101_cudnn7_base .
