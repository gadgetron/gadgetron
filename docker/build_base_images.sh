#!/bin/bash

base_dir=$(pwd)

cd ${base_dir}/base/ubuntu_2004_basic
docker build --network=host -t gadgetron/ubuntu_2004_base .

cd ${base_dir}/base/ubuntu_2004_cuda11_cudnn8
docker build --network=host -t gadgetron/ubuntu_2004_cuda11_cudnn8_base .