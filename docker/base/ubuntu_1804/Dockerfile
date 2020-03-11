FROM ubuntu:18.04

ARG BART_URL=https://github.com/mrirecon/bart
ARG BART_BRANCH=master

RUN apt-get update --quiet && \
    apt-get install --no-install-recommends --no-install-suggests --yes  \
    software-properties-common apt-utils wget build-essential cython3  python3-dev python3-pip libhdf5-serial-dev cmake git-core libboost-all-dev libfftw3-dev h5utils jq hdf\
5-tools liblapack-dev libatlas-base-dev libxml2-dev libfreetype6-dev pkg-config libxslt-dev libarmadillo-dev libace-dev gcc-multilib libgtest-dev liblapacke-dev libplplot-dev libdcmtk-dev sup\
ervisor supervisor net-tools cpio libpugixml-dev jove libopenblas-base libopenblas-dev libpugixml-dev

RUN pip3 install -U pip setuptools
RUN pip3 install numpy scipy Cython tk-tools matplotlib scikit-image opencv_python pydicom scikit-learn h5py
RUN apt-get install --no-install-recommends --no-install-suggests --yes python3-psutil python3-pyxb python3-lxml python3-pil python3-h5py
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends --no-install-suggests --yes python3-tk

RUN pip3 install torch==1.4.0+cpu torchvision==0.5.0+cpu -f https://download.pytorch.org/whl/torch_stable.html
RUN pip3 install --upgrade tensorflow
RUN pip3 install tensorboardx visdom

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
    mkdir build && \
    cd build && \
    make -j $(nproc) && \
    make install

# Ceres
RUN apt-get install --yes libgoogle-glog-dev libeigen3-dev libsuitesparse-dev
RUN cd /opt && \
    wget http://ceres-solver.org/ceres-solver-1.14.0.tar.gz && \
    tar zxf ceres-solver-1.14.0.tar.gz && \
    mkdir ceres-bin && \
    cd ceres-bin && \
    cmake ../ceres-solver-1.14.0 && \
    make -j16 && \
    make install

# since cmake has problems to find python3 and boost-python3
# RUN ln -s /usr/lib/x86_64-linux-gnu/libboost_python-py36.so /usr/lib/x86_64-linux-gnu/libboost_python3.so

# fix the  qhull reentrant problem
# RUN pip uninstall -y scipy

# for embedded python plot, we need agg backend
RUN mkdir -p /root/.config/matplotlib && touch /root/.config/matplotlib/matplotlibrc 
RUN echo "backend : agg" >> /root/.config/matplotlib/matplotlibrc

#Set more environment variables in preparation for Gadgetron installation
ENV GADGETRON_HOME=/usr/local \
    ISMRMRD_HOME=/usr/local

ENV PATH=$PATH:$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib

# Clean up packages.
#RUN  apt-get clean && \
#   rm -rf /var/lib/apt/lists/*
