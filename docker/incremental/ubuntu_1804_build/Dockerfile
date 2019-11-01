FROM gadgetron/ubuntu_1804_cuda100_cudnn7_base

ARG GADGETRON_URL=https://github.com/gadgetron/gadgetron
ARG GADGETRON_BRANCH=master

#ISMRMRD
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd.git && \
    cd ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j $(nproc) && \
    make install

# Install Python interface.
RUN pip3 install gadgetron

#SIEMENS_TO_ISMRMRD
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/siemens_to_ismrmrd.git && \
    cd siemens_to_ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j $(nproc) && \
    make install

#PHILIPS_TO_ISMRMRD
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/philips_to_ismrmrd.git && \
    cd philips_to_ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j $(nproc) && \
    make install

# Clean up packages.
RUN  apt-get clean