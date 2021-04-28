#!/bin/bash

if [ "$EUID" -ne 0 ]
  then echo "Please run as root"
  exit
fi

if [ -z "$(cat /etc/lsb-release | grep "Ubuntu 20.04")" ] && [ -z "$(cat /etc/lsb-release | grep "Ubuntu 18.04")" ]; then
    echo "Error: This install script is intended for Ubuntu 18.04 and 20.04 only"
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
    libxml2-dev \
    libxslt-dev \
    net-tools \
    ninja-build \
    pkg-config \
    python3-dev \
    python3-pip \
    software-properties-common \
    supervisor \
    wget \
    locales

locale-gen en_US.UTF-8
update-locale LANG=en_US.UTF-8

if [ -z "$(cat /etc/lsb-release | grep "Ubuntu 18.04")" ]; then
    # This is NOT ubuntu 18.04, i.e. it is 20.04
    DEBIAN_FRONTEND=noninteractive apt install --no-install-recommends --no-install-suggests --yes \
        googletest \
        googletest-tools \
        librange-v3-dev
else
    # Let's get GCC/G++9
    add-apt-repository --yes --update ppa:ubuntu-toolchain-r/test
    DEBIAN_FRONTEND=noninteractive apt install --no-install-recommends --no-install-suggests --yes \
        g++-9 \
        gcc-9 \
        git

    # Set v9 with higher priority
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 900 --slave /usr/bin/g++ g++ /usr/bin/g++-9

    # install range v3
    mkdir -p /opt/code
    cd /opt/code && \
    git clone https://github.com/ericniebler/range-v3.git --depth 1 --branch 0.11.0 && \
    cd range-v3 && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release ../ -GNinja && \
    ninja && ninja install && cd /opt/code && rm -rf /opt/code/range-v3
    
    # Install Google Test     
    cd /opt/code && \
        git clone https://github.com/google/googletest.git && \
        cd googletest && \
        mkdir build && \
        cd build && \
        cmake -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release ../ && \
        make -j $(nproc) && make install && cd /opt/code && rm -rf /opt/code/googletest
fi

# Install ZFP
mkdir -p /opt/code
    cd /opt/code && \
    git clone https://github.com/LLNL/zfp.git && \
    cd zfp && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    cmake --build . --config Release --parallel $(nproc) && cmake --install .  && cd /opt/code && rm -rf /opt/code/zfp

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
    tk-tools
pip3 install git+https://github.com/ismrmrd/ismrmrd-python.git
pip3 install git+https://github.com/gadgetron/gadgetron-python.git

# If this is an image with CUDA...
if [ -f /usr/local/cuda/bin/nvcc ]; then
    DEBIAN_FRONTEND=noninteractive apt install --no-install-recommends --no-install-suggests --yes libcudnn8-dev
    pip3 install torch==1.7.0+cu110 torchvision==0.8.1+cu110 -f https://download.pytorch.org/whl/torch_stable.html
else
    pip3 install torch==1.7.0+cpu torchvision==0.8.1+cpu -f https://download.pytorch.org/whl/torch_stable.html
fi
