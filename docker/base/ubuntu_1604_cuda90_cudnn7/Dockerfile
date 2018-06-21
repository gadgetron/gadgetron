FROM nvidia/cuda:9.0-cudnn7-devel-ubuntu16.04

RUN apt-get update --quiet && \
    apt-get install --no-install-recommends --no-install-suggests --yes  \
    apt-utils wget build-essential emacs python-dev python-pip python3-dev python3-pip libhdf5-serial-dev cmake git-core \
    libboost-all-dev libfftw3-dev h5utils jq libzmq-dev \
    hdf5-tools liblapack-dev libxml2-dev libfreetype6-dev pkg-config \
    libxslt-dev libarmadillo-dev libace-dev gcc-multilib  \
    libgtest-dev python-dev liblapack-dev liblapacke-dev libopenblas-base libopenblas-dev libplplot-dev libdcmtk-dev \
    libceres-dev supervisor net-tools cpip unzip

#Update gcc to something remotely recent
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get update --quiet
RUN apt-get install gcc-snapshot -y
RUN apt-get -y install gcc-6 g++-6 gfortran
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-6 100 --slave /usr/bin/g++ g++ /usr/bin/g++-6

#Python stuff
RUN pip install --upgrade pip
RUN pip install setuptools
RUN pip install Cython
RUN pip install numpy 
RUN pip install pyxb psutil lxml pydicom
RUN pip install 'h5py==2.5.0' 

RUN apt-get update && apt-get install -y libgtk2.0-dev

RUN pip3 install --upgrade pip
RUN pip3 install -U pip setuptools
RUN apt-get install --no-install-recommends --no-install-suggests --yes python3-psutil python3-pyxb python3-lxml python3-numpy
RUN apt-get install --no-install-recommends --no-install-suggests --yes python3-pil
RUN apt-get install --no-install-recommends --no-install-suggests --yes python3-scipy
RUN apt-get install --no-install-recommends --no-install-suggests --yes python3-configargparse
RUN pip3 install Cython tk-tools matplotlib scikit-image opencv_python pydicom scikit-learn
RUN pip3 uninstall h5py
RUN apt-get install -y python3-h5py

RUN pip3 install --upgrade tensorflow tensorflow-gpu
RUN pip3 install http://download.pytorch.org/whl/cu90/torch-0.4.0-cp35-cp35m-linux_x86_64.whl  
RUN pip3 install torchvision 
RUN pip3 install tensorboardx visdom

# since cmake has problems to find python3 and boost-python3
RUN ln -s /usr/lib/x86_64-linux-gnu/libboost_python-py35.so /usr/lib/x86_64-linux-gnu/libboost_python3.so

# fix the  qhull reentrant problem
# RUN pip uninstall -y scipy

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

ENV PATH=$PATH:/usr/local/cuda-9.0/bin;$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-9.0/lib64:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib

ENV LIBRARY_PATH /usr/local/cuda/lib64/stubs:${LIBRARY_PATH}

# Clean up packages.
#RUN  apt-get clean && \
#   rm -rf /var/lib/apt/lists/*
