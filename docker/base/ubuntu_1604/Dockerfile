FROM ubuntu:16.04

RUN apt-get update --quiet && \
    apt-get install --no-install-recommends --no-install-suggests --yes  \
    wget build-essential emacs python-pip libhdf5-serial-dev cmake git-core \
    libboost-all-dev libfftw3-dev h5utils jq libzmq-dev \
    hdf5-tools liblapack-dev libxml2-dev libfreetype6-dev pkg-config \
    libxslt-dev libarmadillo-dev libace-dev gcc-multilib  \
    libgtest-dev python-dev liblapack-dev liblapacke-dev libatlas-base-dev libatlas-dev libplplot-dev libdcmtk-dev \
    supervisor 

#Python stuff
RUN pip install --upgrade pip
RUN pip install setuptools
RUN pip install Cython
RUN pip install numpy 
RUN pip install pyxb psutil lxml pydicom
RUN pip install 'h5py==2.5.0' 

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
    make install

#BART
RUN cd /opt && \
    wget https://github.com/mrirecon/bart/archive/v0.3.01.tar.gz && \
    tar -xzf v0.3.01.tar.gz && \
    cd bart-0.3.01 && \
    make && \
    ln -s /opt/bart-0.3.01/bart /usr/local/bin/bart 

#Set more environment variables in preparation for Gadgetron installation
ENV GADGETRON_HOME=/usr/local \
    ISMRMRD_HOME=/usr/local

ENV PATH=$PATH:$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib

# Clean up packages.
RUN  apt-get clean && \
   rm -rf /var/lib/apt/lists/*

