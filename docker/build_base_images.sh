#!/bin/bash

base_dir=$(pwd)

cd ${base_dir}/base/ubuntu_1604
docker build -t gadgetron/ubuntu1604_base .

cd ${base_dir}/base/ubuntu_1604_cuda55
docker build -t gadgetron/ubuntu1604_cuda55_base .

cd ${base_dir}/base/ubuntu_1604_cuda75
docker build -t gadgetron/ubuntu1604_cuda75_base .

cd ${base_dir}/base/ubuntu_1604_cuda80
docker build -t gadgetron/ubuntu1604_cuda80_base .

cd ${base_dir}/base/ubuntu_1604_cuda80_cudnn6
docker build -t gadgetron/ubuntu1604_cuda80_cudnn6_base .

cd ${base_dir}/base/ubuntu_1604_cuda80_cudnn7
docker build -t gadgetron/ubuntu1604_cuda80_cudnn7_base .
