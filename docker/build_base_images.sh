#!/bin/bash

base_dir=$(pwd)

cd ${base_dir}/base/ubuntu_1604
docker build --network=host -t gadgetron/ubuntu_1604_base .

cd ${base_dir}/base/ubuntu_1604_cuda55
docker build --network=host -t gadgetron/ubuntu_1604_cuda55_base .

cd ${base_dir}/base/ubuntu_1604_cuda75
docker build --network=host -t gadgetron/ubuntu_1604_cuda75_base .

cd ${base_dir}/base/ubuntu_1604_cuda80
docker build --network=host -t gadgetron/ubuntu_1604_cuda80_base .

cd ${base_dir}/base/ubuntu_1604_cuda80_cudnn6
docker build --network=host -t gadgetron/ubuntu_1604_cuda80_cudnn6_base .

cd ${base_dir}/base/ubuntu_1604_cuda80_cudnn7
docker build --network=host -t gadgetron/ubuntu_1604_cuda80_cudnn7_base .

cd ${base_dir}/base/ubuntu_1604_cuda90_cudnn7
docker build --network=host -t gadgetron/ubuntu_1604_cuda90_cudnn7_base .

cd ${base_dir}/base/ubuntu_1604_cuda92_cudnn7
docker build --network=host -t gadgetron/ubuntu_1604_cuda92_cudnn7_base .

cd ${base_dir}/base/ubuntu_1804
docker build --network=host -t gadgetron/ubuntu_1804_base .

cd ${base_dir}/base/ubuntu_1804_cuda90
docker build --network=host -t gadgetron/ubuntu_1804_cuda90_cudnn7_base .

cd ${base_dir}/base/ubuntu_1804_cuda92_cudnn7
docker build --network=host -t gadgetron/ubuntu_1804_cuda92_cudnn7_base .
