#!/bin/bash

if [ "$EUID" -ne 0 ]
  then echo "Please run as root"
  exit
fi

if [ -z "$(cat /etc/lsb-release | grep "Ubuntu 20.04")" ]; then
    echo "Error: This install script is intended for Ubuntu 20.04 only"
    exit 1
fi

apt update --quiet
DEBIAN_FRONTEND=noninteractive apt install --no-install-recommends --no-install-suggests --yes \
    apt-utils \
    build-essential \
    cmake \
    cpio \
    cython3 \
    gcc-multilib \
    git-core \
    googletest \
    googletest-tools \
    h5utils \
    hdf5-tools \
    jove \
    jq \
    libace-dev \
    libarmadillo-dev \
    libatlas-base-dev \
    libboost-all-dev \
    libcrypto++-dev \
    libfftw3-dev \
    libfreetype6-dev \
    libgtest-dev \
    libhdf5-serial-dev \
    liblapack-dev \
    liblapacke-dev \
    libopenblas-base \
    libopenblas-dev \
    libplplot-dev \
    libpugixml-dev \
    librange-v3-dev \
    libxml2-dev \
    libxslt-dev \
    net-tools \
    ninja-build \
    pkg-config \
    python3-dev \
    python3-pip \
    software-properties-common \
    supervisor \
    wget

pip3 install -U pip setuptools testresources
DEBIAN_FRONTEND=noninteractive apt install --no-install-recommends --no-install-suggests --yes python3-tk

# h5py needs to be recompiled to compile agains HDF5 1.10, which is what we install on Ubuntu 20.04
pip3 install --no-binary=h5py h5py

# Rest of the Python "stuff"
pip3 install numpy scipy Cython tk-tools matplotlib scikit-image opencv_python pydicom scikit-learn sympy Pillow pyxb

# If this is an image with CUDA...
if [ -f /usr/local/cuda/bin/nvcc ]; then
    DEBIAN_FRONTEND=noninteractive apt install --no-install-recommends --no-install-suggests --yes libcudnn8-dev
    pip3 install torch==1.7.0+cu110 torchvision==0.8.1+cu110 -f https://download.pytorch.org/whl/torch_stable.html
else
    pip3 install torch==1.7.0+cpu torchvision==0.8.1+cpu -f https://download.pytorch.org/whl/torch_stable.html
fi
