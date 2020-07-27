FROM nvidia/cuda:11.0-devel-ubuntu20.04

RUN rm -rf /var/lib/apt/lists/*

RUN apt update --quiet && \
    DEBIAN_FRONTEND=noninteractive apt install --no-install-recommends --no-install-suggests --yes  \
    apt-utils software-properties-common wget build-essential cython3 libcrypto++-dev  python3-dev python3-pip libhdf5-serial-dev cmake git-core libboost-all-dev libfftw3-dev h5utils jq hdf5-tools liblapack-dev libatlas-base-dev libxml2-dev libfreetype6-dev pkg-config libxslt-dev libarmadillo-dev libace-dev libgtest-dev liblapacke-dev libplplot-dev gcc-multilib supervisor net-tools cpio libpugixml-dev jove libopenblas-base libopenblas-dev

RUN pip3 install -U pip setuptools testresources
RUN DEBIAN_FRONTEND=noninteractive apt install --no-install-recommends --no-install-suggests --yes python3-tk
RUN pip3 install numpy scipy pandas Cython tk-tools matplotlib seaborn scikit-image opencv_python pydicom scikit-learn sympy Pillow==7.1.2 h5py pyxb
RUN pip3 install torch torchvision
RUN pip3 install onnx onnxruntime

# for embedded python plot, we need agg backend
RUN mkdir -p /root/.config/matplotlib && touch /root/.config/matplotlib/matplotlibrc && echo "backend : agg" >> /root/.config/matplotlib/matplotlibrc

# compile dcmtk
RUN mkdir /opt/code && \
    cd /opt/code && \
    wget https://dicom.offis.de/download/dcmtk/dcmtk365/dcmtk-3.6.5.tar.gz && \
    tar xvf /opt/code/dcmtk-3.6.5.tar.gz && \
    mkdir /opt/code/dcmtk-3.6.5/build && \
    cd /opt/code/dcmtk-3.6.5/build && \
    cmake -DBUILD_SHARED_LIBS=ON ../ && make -j16 && make install && cd /opt/code && \
    rm -rf /opt/code/dcmtk-3.6.5 && rm -rf /opt/code/dcmtk-3.6.5.tar.gz

RUN cd /opt/code && \
    git clone https://github.com/google/googletest.git && \
    cd googletest && \
    mkdir build && \
    cd build && \
    cmake -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release ../ && \
    make -j $(nproc) && make install && cd /opt/code && rm -rf /opt/code/googletest

#Set more environment variables in preparation for Gadgetron installation
ENV GADGETRON_HOME=/usr/local \
    ISMRMRD_HOME=/usr/local

ENV PATH=$PATH:/usr/local/cuda/bin;$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib

ENV LIBRARY_PATH /usr/local/cuda/lib64/stubs:${LIBRARY_PATH}
