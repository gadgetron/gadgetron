FROM gadgetron/ubuntu_1804_cuda100_cudnn7_base

ARG GADGETRON_URL=https://github.com/gadgetron/gadgetron
ARG GADGETRON_BRANCH=master

# RUN apt-get -y install gcc-6 g++-6 gfortran

#ISMRMRD
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd.git && \
    cd ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j $(nproc) && \
    make install

#Fix compiler error with CUDA for now
RUN sed -i "s;#error -- unsupported GNU version! gcc versions later than 4.9 are not supported!;//#error -- unsupported GNU version! gcc versions later than 4.9 are not supported!;g" /usr/local/cuda/include/host_config.h

#GADGETRON
RUN cd /opt/code && \
    git clone ${GADGETRON_URL} --branch ${GADGETRON_BRANCH} --single-branch && \
    cd gadgetron && \
    mkdir build && \
    cd build && \
    cmake -DBUILD_WITH_PYTHON3=ON ../ && \
    make -j $(nproc) && \
    make install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.gadgetron.sha1 --value `git rev-parse HEAD` && \
    cp ${GADGETRON_HOME}/share/gadgetron/config/gadgetron.xml.example ${GADGETRON_HOME}/share/gadgetron/config/gadgetron.xml && \
    cp /opt/code/gadgetron/docker/start_supervisor /opt/ && \
    cp /opt/code/gadgetron/docker/supervisord.conf /opt/

#HASH for ISMRMRD
RUN cd /opt/code/ismrmrd && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.ismrmrd.sha1 --value `git rev-parse HEAD` 

#ISMRMRD PYTHON API
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd-python.git &&  \
    cd ismrmrd-python && \
    python3 setup.py install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.ismrmrdpython.sha1 --value `git rev-parse HEAD` 


#ISMRMRD PYTHON TOOLS
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd-python-tools.git &&  \
    cd ismrmrd-python-tools && \
    python3 setup.py install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.ismrmrdpythontools.sha1 --value `git rev-parse HEAD` 

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
