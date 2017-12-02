FROM nvidia/cuda:8.0-cudnn7-devel-ubuntu16.04

RUN apt-get update --quiet && \
    apt-get install --no-install-recommends --no-install-suggests --yes  \
    apt-utils wget build-essential emacs python-pip libhdf5-serial-dev cmake git-core \
    libboost-all-dev libfftw3-dev h5utils jq libzmq-dev \
    hdf5-tools liblapack-dev libxml2-dev libfreetype6-dev pkg-config \
    libxslt-dev libarmadillo-dev libace-dev gcc-multilib  \
    libgtest-dev python-dev liblapack-dev liblapacke-dev libatlas-base-dev libatlas-dev libplplot-dev libdcmtk-dev \
    supervisor net-tools

#Python stuff
RUN pip install --upgrade pip
RUN pip install setuptools
RUN pip install Cython
RUN pip install numpy 
RUN pip install pyxb psutil lxml pydicom
RUN pip install 'h5py==2.5.0' 

RUN apt-get update && apt-get install -y python python-setuptools python-pip python-dev build-essential libgtk2.0-dev
RUN easy_install -U pip
RUN apt install --yes python-scipy
RUN pip install tensorflow-gpu matplotlib opencv-python scikit-image

# fix the  qhull reentrant problem
# RUN pip uninstall -y scipy

#OpenBLAS with OpenMP
RUN cd /opt && \
    mkdir debsource && \
    cd debsource && \
    apt-get --no-install-recommends --no-install-suggests --yes build-dep libopenblas-base && \
    apt-get install --no-install-recommends --no-install-suggests --yes build-essential fakeroot devscripts && \
    apt-get source libopenblas-base && \
    cd openblas-0.2.18/ && \
    sed -i "s/NO_WARMUP=1/NO_WARMUP=1 OPENMP=1/g" debian/rules && \
    debchange -i "Compiling with OpenMP support" && \
    debuild -us -uc -i -I && \
    debi && \
    update-alternatives --set libblas.so.3 /usr/lib/openblas-base/libblas.so.3 && \
    cd /opt && \
    rm -rf debsource

#ZFP
RUN cd /opt && \
    git clone https://github.com/hansenms/ZFP.git && \
    cd ZFP && \
    mkdir lib && \
    make && \
    make shared && \
    make -j 16 install

#BART
RUN cd /opt && \
    wget https://github.com/mrirecon/bart/archive/v0.3.01.tar.gz && \
    tar -xzf v0.3.01.tar.gz && \
    cd bart-0.3.01 && \
    make -j 16 && \
    ln -s /opt/bart-0.3.01/bart /usr/local/bin/bart 

#Set more environment variables in preparation for Gadgetron installation
ENV GADGETRON_HOME=/usr/local \
    ISMRMRD_HOME=/usr/local

ENV PATH=$PATH:$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib

ENV LIBRARY_PATH /usr/local/cuda/lib64/stubs:${LIBRARY_PATH}

# Clean up packages.
#RUN  apt-get clean && \
#   rm -rf /var/lib/apt/lists/*
