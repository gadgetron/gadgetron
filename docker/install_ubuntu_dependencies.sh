#!/bin/bash

if [ "$EUID" -ne 0 ]; then
  echo "Please run as root"
  exit
fi

if [ -z "$(cat /etc/lsb-release | grep "Ubuntu 20.04")" ] ; then
  echo "Error: This install script is intended for Ubuntu  20.04 only"
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
  h5utils \
  hdf5-tools \
  jove \
  jq \
  libace-dev \
  libarmadillo-dev \
  libatlas-base-dev \
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
  libxml2-dev \
  libxslt-dev \
  librocksdb-dev \
  net-tools \
  ninja-build \
  pkg-config \
  python3-dev \
  python3-pip \
  software-properties-common \
  supervisor \
  wget \
  googletest \
  googletest-tools \
  librange-v3-dev \
  nlohmann-json3-dev \
  libboost-all-dev

# Install ZFP
mkdir -p /opt/code
cd /opt/code &&
  git -c advice.detachedHead=false clone --branch 0.5.5 --single-branch https://github.com/LLNL/zfp.git \
    cd zfp &&
  mkdir build &&
  cd build &&
  cmake ../ &&
  cmake --build . --config Release --parallel $(nproc) && cmake --install . && cd /opt/code && rm -rf /opt/code/zfp

pip3 install -U pip setuptools testresources
DEBIAN_FRONTEND=noninteractive apt install --no-install-recommends --no-install-suggests --yes python3-tk

# h5py needs to be recompiled to compile agains HDF5 1.10, which is what we install on Ubuntu 20.04
pip3 install --no-binary=h5py h5py

# Rest of the Python "stuff"
pip3 install \
  Cython \
  matplotlib \
  numpy \
  opencv_python \
  pydicom \
  Pillow \
  pyxb \
  scikit-image \
  scikit-learn \
  scipy \
  sympy \
  tk-tools \
  junitparser

env LC_ALL=C.UTF-8 LANG=C.UTF-8 pip3 install git+https://github.com/ismrmrd/ismrmrd-python.git

pip3 install git+https://github.com/gadgetron/gadgetron-python.git

wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.4-linux-x86_64.tar.gz
tar zxvf julia-1.6.4-linux-x86_64.tar.gz -C /opt/
cp /opt/julia-1.6.4/* /usr/local
rm julia-1.6.4-linux-x86_64.tar.gz
rm -rf /opt/julia-1.6.4

julia -e "Import Pkg; Pkg.add(\"https://github.com/gadgetron/Gadgetron.jl.git\");
        Pkg.add(\"https://github.com/gadgetron/GadgetronExamples.jl.git\");"

# If this is an image with CUDA...
if [ -f /usr/local/cuda/bin/nvcc ]; then
  DEBIAN_FRONTEND=noninteractive apt install --no-install-recommends --no-install-suggests --yes libcudnn8-dev
  pip3 install torch==1.7.0+cu110 torchvision==0.8.1+cu110 -f https://download.pytorch.org/whl/torch_stable.html
else
  pip3 install torch==1.7.0+cpu torchvision==0.8.1+cpu -f https://download.pytorch.org/whl/torch_stable.html
fi
