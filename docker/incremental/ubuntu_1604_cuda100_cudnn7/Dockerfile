FROM gadgetron/ubuntu_1604_cuda100_cudnn7_base

ARG GADGETRON_URL=https://github.com/gadgetron/gadgetron
ARG GADGETRON_BRANCH=master

#ISMRMRD
RUN mkdir /opt/code
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd.git && \
    cd ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j $(nproc) && \
    make install

#GADGETRON
RUN cd /opt/code && \
    git clone ${GADGETRON_URL} --branch ${GADGETRON_BRANCH} --single-branch && \
    cd gadgetron && \
    mkdir build && \
    cd build && \
    cmake -DCUDA_NVCC_FLAGS="-ccbin gcc-7" ../ && \
    make -j $(nproc) && \
    make install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.gadgetron.sha1 --value `git rev-parse HEAD` && \
    cp /opt/code/gadgetron/docker/start_supervisor /opt/ && \
    cp /opt/code/gadgetron/docker/supervisord.conf /opt/

#HASH for ISMRMRD
RUN cd /opt/code/ismrmrd && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.ismrmrd.sha1 --value `git rev-parse HEAD` 

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
    make install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.siemens_to_ismrmrd.sha1 --value `git rev-parse HEAD`

#PHILIPS_TO_ISMRMRD
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/philips_to_ismrmrd.git && \
    cd philips_to_ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j $(nproc) && \
    make install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.philips_to_ismrmrd.sha1 --value `git rev-parse HEAD` 

# Clean up packages.
RUN  apt-get clean && \
   rm -rf /var/lib/apt/lists/*

CMD ["/opt/start_supervisor"]
