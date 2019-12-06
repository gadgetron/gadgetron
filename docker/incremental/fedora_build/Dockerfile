FROM fedora:30

RUN dnf update -y
RUN dnf install -y \
    git-core wget make cmake gcc gcc-c++ \
    boost-devel boost-python3-devel armadillo-devel libxml2-devel libxslt-devel pugixml-devel fftw-devel \
    python3-devel python3-h5py python3-scipy openblas-openmp python3-pip
RUN pip3 install PyXB

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

# Install Python interface.
RUN pip3 install gadgetron

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

RUN mkdir /var/run/sshd
RUN echo 'root:root' | chpasswd
RUN sed -i 's/PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config

# SSH login fix. Otherwise user is kicked off after login
RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd

ENV NOTVISIBLE "in users profile"
RUN echo "export VISIBLE=now" >> /etc/profile

# 22 for ssh server. 7777 for gdb server.
EXPOSE 22 7777

RUN useradd -ms /bin/bash debugger
RUN echo 'debugger:pwd' | chpasswd

########################################################
# Add custom packages and development environment here
########################################################

########################################################

CMD ["/usr/sbin/sshd", "-D"]

