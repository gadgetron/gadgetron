FROM ubuntu:18.04

ARG BART_URL=https://github.com/mrirecon/bart
ARG BART_BRANCH=master

RUN apt-get update --quiet && \
    apt-get install --no-install-recommends --no-install-suggests --yes  \
    software-properties-common apt-utils wget build-essential cython emacs python-dev python-pip python3-dev python3-pip libhdf5-serial-dev cmake git-core libboost-all-dev libfftw3-dev h5utils jq hdf\
5-tools liblapack-dev libxml2-dev libfreetype6-dev pkg-config libxslt-dev libarmadillo-dev libace-dev gcc-multilib libgtest-dev python3-dev liblapack-dev liblapacke-dev libplplot-dev libdcmtk-dev sup\
ervisor cmake-curses-gui neofetch supervisor net-tools cpio gpg-agent

# install cuda 9
RUN apt-get install --no-install-recommends --no-install-suggests --yes linux-image-extra-virtual
RUN apt-get install --no-install-recommends --no-install-suggests --yes linux-source
# RUN apt-get source linux-image-$(uname -r)
RUN apt-get install --no-install-recommends --no-install-suggests --yes linux-headers-$(uname -r)

# RUN add-apt-repository ppa:graphics-drivers/ppa
# RUN apt update --quiet && apt upgrade -y
# RUN DEBIAN_FRONTEND=noninteractive apt install --no-install-recommends --no-install-suggests --yes nvidia-390

RUN mkdir /opt/code
RUN cd /opt/code

RUN wget https://developer.nvidia.com/compute/cuda/9.0/Prod/local_installers/cuda-repo-ubuntu1704-9-0-local_9.0.176-1_amd64-deb
RUN dpkg -i cuda-repo-ubuntu1704-9-0-local_9.0.176-1_amd64-deb
RUN apt-key add /var/cuda-repo-9-0-local/7fa2af80.pub
RUN apt-get update --quiet 
RUN DEBIAN_FRONTEND=noninteractive apt-get install cuda --no-install-recommends --no-install-suggests --yes

RUN wget https://developer.nvidia.com/compute/cuda/9.0/Prod/patches/1/cuda-repo-ubuntu1704-9-0-local-cublas-performance-update_1.0-1_amd64-deb
RUN dpkg -i cuda-repo-ubuntu1704-9-0-local-cublas-performance-update_1.0-1_amd64-deb
RUN apt-get update --quiet 
RUN DEBIAN_FRONTEND=noninteractive apt-get upgrade --yes cuda

RUN wget https://developer.nvidia.com/compute/cuda/9.0/Prod/patches/2/cuda-repo-ubuntu1704-9-0-local-cublas-performance-update-2_1.0-1_amd64-deb
RUN dpkg -i cuda-repo-ubuntu1704-9-0-local-cublas-performance-update-2_1.0-1_amd64-deb
RUN apt-get update --quiet 
RUN DEBIAN_FRONTEND=noninteractive apt-get upgrade --yes cuda

RUN wget https://developer.nvidia.com/compute/cuda/9.0/Prod/patches/3/cuda-repo-ubuntu1704-9-0-local-cublas-performance-update-3_1.0-1_amd64-deb
RUN dpkg -i cuda-repo-ubuntu1704-9-0-local-cublas-performance-update-3_1.0-1_amd64-deb
RUN apt-get update --quiet 
RUN DEBIAN_FRONTEND=noninteractive apt-get upgrade --yes cuda

COPY /docker/cudnn-9.0-linux-x64-v7.tgz /opt/code/

RUN cd /opt/code && \
    tar -xvf ./cudnn-9.0-linux-x64-v7.tgz && \
    cp cuda/include/cudnn.h /usr/local/cuda/include && \
    cp cuda/lib64/libcudnn* /usr/local/cuda/lib64 && \
    chmod a+r /usr/local/cuda/include/cudnn.h /usr/local/cuda/lib64/libcudnn*

RUN pip3 install --upgrade pip
RUN pip3 install -U pip setuptools
RUN apt-get install --no-install-recommends --no-install-suggests --yes python3-psutil python3-pyxb python3-lxml
RUN apt-get install --no-install-recommends --no-install-suggests --yes python3-pil
RUN apt-get install --no-install-recommends --no-install-suggests --yes python3-scipy
RUN apt-get install --no-install-recommends --no-install-suggests --yes python3-configargparse
RUN pip3 install numpy==1.15.4 Cython tk-tools matplotlib scikit-image opencv_python pydicom scikit-learn
RUN pip3 uninstall h5py
RUN apt-get install -y python3-h5py

RUN pip3 install --upgrade tensorflow tensorflow-gpu
RUN pip3 install torch
RUN pip3 install torchvision 
RUN pip3 install tensorboardx visdom

# since cmake has problems to find python3 and boost-python3
# RUN ln -s /usr/lib/x86_64-linux-gnu/libboost_python-py36.so /usr/lib/x86_64-linux-gnu/libboost_python3.so

# fix the  qhull reentrant problem
# RUN pip uninstall -y scipy

#OpenBLAS with OpenMP
RUN cd /opt && \
    mkdir debsource && \
    cd debsource && \
    apt-get --no-install-recommends --no-install-suggests --yes build-dep libopenblas-base && \
    apt-get install --no-install-recommends --no-install-suggests --yes build-essential fakeroot devscripts && \
    apt-get source libopenblas-base && \
    cd openblas-0.2.20+ds/ && \
    sed -i "s/NO_WARMUP=1/NO_WARMUP=1 OPENMP=1/g" debian/rules && \
    debchange -i "Compiling with OpenMP support" && \
    debuild -us -uc -i -I && \
    debi && \
    cd /opt && \
    rm -rf debsource

#Install Openblas
RUN rm /usr/lib/x86_64-linux-gnu/libblas.so
RUN rm /usr/lib/x86_64-linux-gnu/libblas.so.3
RUN ln -s /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3 /usr/lib/x86_64-linux-gnu/libblas.so
RUN ln -s /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3 /usr/lib/x86_64-linux-gnu/libblas.so.3

RUN rm /usr/lib/x86_64-linux-gnu/liblapack.so
RUN rm /usr/lib/x86_64-linux-gnu/liblapack.so.3
RUN ln -s /usr/lib/x86_64-linux-gnu/openblas/liblapack.so.3 /usr/lib/x86_64-linux-gnu/liblapack.so
RUN ln -s /usr/lib/x86_64-linux-gnu/openblas/liblapack.so.3 /usr/lib/x86_64-linux-gnu/liblapack.so.3

#ZFP
RUN cd /opt && \
    git clone https://github.com/hansenms/ZFP.git && \
    cd ZFP && \
    mkdir lib && \
    make && \
    make shared && \
    make -j $(nproc) install

# BART
RUN cd /opt/code && \
    git clone ${BART_URL} --branch ${BART_BRANCH} --single-branch && \
    cd bart && \
    mkdir build && \
    cd build && \
    cmake .. -DBART_FPIC=ON -DBART_ENABLE_MEM_CFL=ON -DBART_REDEFINE_PRINTF_FOR_TRACE=ON -DBART_LOG_BACKEND=ON -DBART_LOG_GADGETRON_BACKEND=ON && \
    make -j $(nproc) && \
    make install

# ceres
RUN apt-get install --yes libgoogle-glog-dev libeigen3-dev libsuitesparse-dev
RUN cd /opt && \
    wget http://ceres-solver.org/ceres-solver-1.14.0.tar.gz && \
    tar zxf ceres-solver-1.14.0.tar.gz && \
    mkdir ceres-bin && \
    cd ceres-bin && \
    cmake ../ceres-solver-1.14.0 && \
    make -j20 && \
    make install

#Set more environment variables in preparation for Gadgetron installation
ENV GADGETRON_HOME=/usr/local \
    ISMRMRD_HOME=/usr/local

ENV PATH=$PATH:/usr/local/cuda-9.0/bin;$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-9.0/lib64:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib

ENV LIBRARY_PATH /usr/local/cuda/lib64/stubs:${LIBRARY_PATH}

# Clean up packages.
#RUN  apt-get clean && \
#   rm -rf /var/lib/apt/lists/*
