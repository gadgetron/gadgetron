FROM nvidia/cuda:10.0-cudnn7-devel-ubuntu18.04

ARG BART_URL=https://github.com/mrirecon/bart
ARG BART_BRANCH=master

RUN rm -rf /var/lib/apt/lists/*

RUN  sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list
RUN apt-get update --quiet && \
    apt-get install --no-install-recommends --no-install-suggests --yes  \
    software-properties-common apt-utils wget build-essential cython emacs \
    python-dev python-pip python3-dev python3-pip libhdf5-serial-dev cmake \
    git-core libboost-all-dev libfftw3-dev h5utils jq hdf5-tools libxml2-dev \
    libfreetype6-dev pkg-config libxslt-dev gcc-multilib libgtest-dev python3-dev \
    libplplot-dev libdcmtk-dev supervisor cmake-curses-gui neofetch supervisor net-tools \
    cpio gpg-agent libpugixml-dev librange-v3-dev

RUN apt-get install --no-install-recommends --no-install-suggests --yes libopenblas-dev libarmadillo-dev

#Python stuff
RUN apt-get update && apt-get install -y libgtk2.0-dev

RUN pip3 install --upgrade pip
RUN pip3 install -U pip setuptools
RUN pip3 install h5py psutil pyxb lxml Pillow configargparse numpy Cython scipy tk-tools matplotlib scikit-image opencv_python pydicom scikit-learn tqdm

RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends --no-install-suggests --yes python3-tk

RUN pip3 install --upgrade tensorflow
RUN pip3 install torch===1.4.0 torchvision===0.5.0 -f https://download.pytorch.org/whl/torch_stable.html
RUN pip3 install tensorboardx visdom

# for embedded python plot, we need agg backend
RUN mkdir -p /root/.config/matplotlib && touch /root/.config/matplotlib/matplotlibrc 
RUN echo "backend : agg" >> /root/.config/matplotlib/matplotlibrc

RUN mkdir /opt/code

# ZFP
RUN cd /opt && \
    git clone https://github.com/hansenms/ZFP.git && \
    cd ZFP && \
    mkdir lib && \
    make && \
    make shared && \
    make -j $(nproc) install

# BART
RUN apt-get install --no-install-recommends --no-install-suggests --yes liblapacke-dev
RUN cd /opt/code && \
    git clone ${BART_URL} --branch ${BART_BRANCH} --single-branch && \
    cd bart && \
    make -j $(nproc) && \
    make install

#Set more environment variables in preparation for Gadgetron installation
ENV GADGETRON_HOME=/usr/local \
    ISMRMRD_HOME=/usr/local

ENV PATH=$PATH:/usr/local/cuda/bin;$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib

ENV LIBRARY_PATH /usr/local/cuda/lib64/stubs:${LIBRARY_PATH}
