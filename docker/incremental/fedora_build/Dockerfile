FROM fedora:29

RUN dnf update -y
RUN dnf install -y \
    git-core wget make cmake gcc gcc-c++ \
    boost-devel boost-python3-devel armadillo-devel libxml2-devel libxslt-devel pugixml-devel fftw-devel \
    python3-devel python3-h5py python3-scipy python3-PyXB

RUN mkdir /opt/code

#ISMRMRD
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd.git && \
    cd ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j $(nproc) && \
    make install

#ISMRMRD PYTHON API
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd-python.git &&  \
    cd ismrmrd-python && \
    python3 setup.py install

#ISMRMRD PYTHON TOOLS
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd-python-tools.git &&  \
    cd ismrmrd-python-tools && \
    python3 setup.py install

#SIEMENS_TO_ISMRMRD
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/siemens_to_ismrmrd.git && \
    cd siemens_to_ismrmrd && \
    mkdir build && \
    cd build && \
    cmake -DBUILD_DYNAMIC=true ../ && \
    make -j $(nproc) && \
    make install

ENV GADGETRON_HOME=/usr/local \
    ISMRMRD_HOME=/usr/local

ENV LD_LIBRARY_PATH=$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib
