#!/bin/bash

base_dir=$(pwd)

cd ${base_dir}/base/ubuntu_1604_cuda100_cudnn7
docker build --network=host -t gadgetron/ubuntu_1604_cuda100_cudnn7_base .
cd ${base_dir}/base/ubuntu_1804_cuda100_cudnn7
docker build --network=host -t gadgetron/ubuntu_1804_cuda100_cudnn7_base .

cd ${base_dir}
