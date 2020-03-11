FROM nvidia/cuda:10.0-cudnn7-devel-ubuntu16.04

ARG BART_URL=https://github.com/mrirecon/bart
ARG BART_BRANCH=master

#Fix default ceres package on 16.04 is broken

RUN apt-get update --quiet && \
    apt-get install --no-install-recommends --no-install-suggests --yes software-properties-common

RUN apt-get update --quiet && \
    apt-get install --no-install-recommends --no-install-suggests --yes  \
    apt-utils wget build-essential emacs python-dev python-pip libhdf5-serial-dev git-core \
    libfftw3-dev h5utils jq libzmq-dev \
    hdf5-tools liblapack-dev libxml2-dev libfreetype6-dev pkg-config \
    libxslt-dev gcc-multilib  \
    libgtest-dev python-dev liblapack-dev liblapacke-dev libplplot-dev libdcmtk-dev \
    supervisor cpio net-tools libpugixml-dev

RUN apt-get install --no-install-recommends --no-install-suggests --yes \
    libopenblas-base libopenblas-dev

#Update cmake to something not entirely ancient
RUN add-apt-repository ppa:litenstein/cmake310-xenial
RUN apt-get update
RUN apt-get -y install cmake

#Update gcc to something remotely recent
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get update
RUN apt-get -y install gcc-7 g++-7
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 100 --slave /usr/bin/g++ g++ /usr/bin/g++-7

RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get update
RUN apt-get install -y python3.6
RUN apt-get install -y python3.6-dev
RUN apt-get install -y python3.6-venv

RUN cd /tmp
RUN wget https://bootstrap.pypa.io/get-pip.py
RUN python3.6 get-pip.py
RUN rm /usr/bin/python3
RUN ln -s /usr/bin/python3.6 /usr/bin/python3

#Python stuff
RUN pip3 install --upgrade pip
RUN pip3 install -U pip setuptools
RUN pip3 install h5py psutil pyxb lxml Pillow configargparse numpy==1.15.4 Cython scipy tk-tools matplotlib==2.2.3 scikit-image opencv_python pydicom scikit-learn

RUN apt-get install -y python3.6-tk

RUN pip3 install --upgrade tensorflow
RUN pip3 install torch===1.4.0 torchvision===0.5.0 -f https://download.pytorch.org/whl/torch_stable.html
RUN pip3 install tensorboardx visdom tqdm

# for embedded python plot, we need agg backend
RUN mkdir -p /root/.config/matplotlib && touch /root/.config/matplotlib/matplotlibrc 
RUN echo "backend : agg" >> /root/.config/matplotlib/matplotlibrc

# Compile boost
RUN cd /opt && \
    wget https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.gz && \
    tar -xzvf boost_1_69_0.tar.gz && \
    cd boost_1_69_0 && \
    ./bootstrap.sh --with-python=/usr/bin/python3 --prefix=/usr && \
    ./b2 -j$(nproc) address-model=64 stage && \
    ./b2 install

# since cmake has problems to find python3 and boost-python3
RUN ln -s /usr/lib/libboost_python36.so /usr/lib/libboost_python3.so

# Compile armadillo
RUN cd /opt && \
    wget http://sourceforge.net/projects/arma/files/armadillo-9.200.8.tar.xz && \
    tar xvf armadillo-9.200.8.tar.xz && \
    cd armadillo-9.200.8 && \
    mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release ../ && \
    make -j$(nproc) && \
    make PREFIX=/usr/lib/x86_64-linux-gnu install

# Compile PugiXML
RUN cd /opt && \
    git clone https://github.com/zeux/pugixml.git --branch v1.9 && cd pugixml && \
    mkdir build && cd build &&\
    cmake -DCMAKE_BUILD_TYPE=Release .. &&  \
    make -j$(nproc) && \
    make install
    
# fix the  qhull reentrant problem
# RUN pip uninstall -y scipy

#ZFP
RUN cd /opt && \
    git clone https://github.com/hansenms/ZFP.git && \
    cd ZFP && \
    mkdir lib && \
    make && \
    make shared && \
    make -j $(nproc) install

# BART
RUN cd /opt && \
    git clone ${BART_URL} --branch ${BART_BRANCH} --single-branch && \
    cd bart && \
    make -j $(nproc) && \
    make install

RUN apt-get clean

#Set more environment variables in preparation for Gadgetron installation
ENV GADGETRON_HOME=/usr/local \
    ISMRMRD_HOME=/usr/local

ENV PATH=$PATH:/usr/local/cuda/bin;$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib

ENV LIBRARY_PATH /usr/local/cuda/lib64/stubs:${LIBRARY_PATH}

# Clean up packages.
#RUN  apt-get clean && \
#   rm -rf /var/lib/apt/lists/*
