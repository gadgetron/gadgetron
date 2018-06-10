FROM nvidia/cuda:9.2-cudnn7-devel-ubuntu18.04

RUN apt-get update --quiet && \
    apt-get install --no-install-recommends --no-install-suggests --yes  \
    -y apt-utils wget build-essential emacs python-dev python-pip python3-dev python3-pip libhdf5-serial-dev cmake git-core libboost-all-dev libfftw3-dev h5utils jq hdf5-tools liblapack-dev libxml2-dev libfreetype6-dev pkg-config libxslt-dev libarmadillo-dev libace-dev gcc-multilib libgtest-dev python3-dev liblapack-dev liblapacke-dev libplplot-dev libdcmtk-dev supervisor cmake-curses-gui neofetch supervisor net-tools

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
RUN pip3 install Cython
RUN apt-get install -y python3-psutil python3-pyxb python3-lxml
RUN pip3 install pydicom
RUN apt-get install -y python3-h5py
RUN apt-get install -y python3-numpy python3-configargparse
RUN apt install --yes python3-scipy

RUN pip3 install --upgrade tensorflow tensorflow-gpu
RUN pip3 install http://download.pytorch.org/whl/cu91/torch-0.4.0-cp36-cp36m-linux_x86_64.whl 
RUN pip3 install torchvision scipy matplotlib scikit-learn scikit-image jupyter
RUN pip3 install tensorboardx visdom opencv-python

# since cmake has problems to find python3 and boost-python3
RUN ln -s /usr/lib/x86_64-linux-gnu/libboost_python-py36.so /usr/lib/x86_64-linux-gnu/libboost_python3.so

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

ENV PATH=$PATH:/usr/local/cuda-9.2/bin;$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-9.2/lib64:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib

ENV LIBRARY_PATH /usr/local/cuda/lib64/stubs:${LIBRARY_PATH}

# Clean up packages.
#RUN  apt-get clean && \
#   rm -rf /var/lib/apt/lists/*
