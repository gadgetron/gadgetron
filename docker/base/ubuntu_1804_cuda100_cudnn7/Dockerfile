FROM nvidia/cuda:10.0-cudnn7-devel-ubuntu18.04

ARG BART_URL=https://github.com/mrirecon/bart
ARG BART_BRANCH=master

RUN rm -rf /var/lib/apt/lists/*
RUN apt-get update --quiet && \
    apt-get install --no-install-recommends --no-install-suggests --yes  \
    software-properties-common apt-utils wget build-essential cython emacs python-dev python-pip python3-dev python3-pip libhdf5-serial-dev cmake git-core libboost-all-dev libfftw3-dev h5utils jq hdf\
5-tools libatlas-base-dev liblapack-dev libxml2-dev libfreetype6-dev pkg-config libxslt-dev libarmadillo-dev libace-dev gcc-multilib libgtest-dev python3-dev liblapack-dev liblapacke-dev libplplot-dev libdcmtk-dev sup\
ervisor cmake-curses-gui neofetch supervisor net-tools cpio gpg-agent

RUN apt-get install --no-install-recommends --no-install-suggests --yes libopenblas-base libopenblas-dev

#Python stuff
RUN apt-get update && apt-get install -y libgtk2.0-dev

RUN pip3 install --upgrade pip
RUN pip3 install -U pip setuptools
RUN pip3 install h5py psutil pyxb lxml Pillow configargparse numpy==1.15.4 Cython scipy tk-tools matplotlib==2.2.3 scikit-image opencv_python pydicom scikit-learn tqdm
RUN apt-get install --no-install-recommends --no-install-suggests --yes python3-psutil python3-pyxb python3-lxml python3-pil python3-h5py
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends --no-install-suggests --yes python3-tk

RUN pip3 install --upgrade tensorflow
RUN pip3 install https://download.pytorch.org/whl/cu100/torch-1.1.0-cp36-cp36m-linux_x86_64.whl
RUN pip3 install https://download.pytorch.org/whl/cu100/torchvision-0.3.0-cp36-cp36m-linux_x86_64.whl
RUN pip3 install tensorboardx visdom

# since cmake has problems to find python3 and boost-python3
# RUN ln -s /usr/lib/x86_64-linux-gnu/libboost_python-py36.so /usr/lib/x86_64-linux-gnu/libboost_python3.so

# fix the  qhull reentrant problem
# RUN pip uninstall -y scipy

# for embedded python plot, we need agg backend
RUN mkdir -p /root/.config/matplotlib && touch /root/.config/matplotlib/matplotlibrc 
RUN echo "backend : agg" >> /root/.config/matplotlib/matplotlibrc

RUN mkdir /opt/code

#ZFP
RUN cd /opt && \
    git clone https://github.com/hansenms/ZFP.git && \
    cd ZFP && \
    mkdir lib && \
    make && \
    make shared && \
    make -j $(nproc) install

#BART
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

ENV PATH=$PATH:/usr/local/cuda/bin;$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib

ENV LIBRARY_PATH /usr/local/cuda/lib64/stubs:${LIBRARY_PATH}

# Clean up packages.
#RUN  apt-get clean && \
#   rm -rf /var/lib/apt/lists/*
