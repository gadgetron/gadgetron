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

#Fix compiler error with CUDA for now
RUN sed -i "s;#error -- unsupported GNU version! gcc versions later than 4.9 are not supported!;//#error -- unsupported GNU version! gcc versions later than 4.9 are not supported!;g" /usr/local/cuda/include/host_config.h

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

RUN apt-get install clang -y
RUN apt-get install google-perftools libgoogle-perftools-dev -y

# Clean up packages.
RUN  apt-get clean && \
   rm -rf /var/lib/apt/lists/*


#RUN echo 'root:root' | chpasswd
#RUN sed -i 's/#PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config
#
## SSH login fix. Otherwise user is kicked off after login
#RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd
#RUN sed -i 's/#PubkeyAuthentication yes/PubkeyAuthentication yes/' /etc/ssh/sshd_config
#
#ENV NOTVISIBLE "in users profile"
#RUN echo "export VISIBLE=now" >> /etc/profile
#RUN useradd build -p build


#EXPOSE 22 7777
#CMD ["/usr/sbin/sshd", "-D"]