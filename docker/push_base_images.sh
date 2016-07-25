#!/bin/bash

base_dir=$(pwd)

docker push gadgetron/ubuntu1404_base
docker push gadgetron/ubuntu_1404_cuda55_base
docker push gadgetron/ubuntu_1404_cuda75_base
docker push gadgetron/ubuntu1604_base
docker push gadgetron/ubuntu_1604_cuda55_base
docker push gadgetron/ubuntu_1604_cuda75_base

