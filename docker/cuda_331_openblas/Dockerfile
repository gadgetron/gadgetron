FROM ubuntu:14.04

ENV CUDA_DRIVER 331.113
ENV CUDA_DRIVER_INSTALL http://us.download.nvidia.com/XFree86/Linux-x86_64/${CUDA_DRIVER}/NVIDIA-Linux-x86_64-${CUDA_DRIVER}.run

ENV CUDA_DOWNLOAD_PATH http://developer.download.nvidia.com/compute/cuda/5_5/rel/installers 
ENV CUDA_RUN cuda_5.5.22_linux_64.run
ENV CUDA_INSTALL cuda-linux64-rel-5.5.22-16488124.run

# Update and install minimal.
RUN \
  apt-get update \
            --quiet \
  && apt-get install \
            --yes \
            --no-install-recommends \
            --no-install-suggests \
       build-essential \
       module-init-tools \
       wget

# Install CUDA DRIVER
RUN wget \
      $CUDA_DRIVER_INSTALL -P /tmp --no-verbose \
      && chmod +x /tmp/NVIDIA-Linux-x86_64-${CUDA_DRIVER}.run \
      && /tmp/NVIDIA-Linux-x86_64-${CUDA_DRIVER}.run -s -N --no-kernel-module \
      && rm -rf /tmp/*

#UNPACK THE CUDA INSTALLER
RUN cd /opt && \
  wget ${CUDA_DOWNLOAD_PATH}/${CUDA_RUN} --no-verbose && \
  chmod +x *.run && \
  mkdir nvidia_installers && \
  ./${CUDA_RUN} -extract=`pwd`/nvidia_installers && \
  rm ${CUDA_RUN}

#INSTALL CUDA LIBRARIES, etc. 
RUN cd /opt/nvidia_installers && \
    ./${CUDA_INSTALL} -noprompt && \
    cd ../ && \
    rm -rf nvidia_installers

#SET SOME ENV variables for CUDA
ENV CUDA_PATH=/usr/local/cuda \
    LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64 \
    PATH=$PATH:/usr/local/cuda/bin

RUN apt-get install --no-install-recommends --no-install-suggests --yes  \
    emacs python-pip cython libhdf5-serial-dev cmake git-core \
    libboost-all-dev libfftw3-dev h5utils \
    hdf5-tools python-dev python-numpy liblapack-dev libxml2-dev \
    libxslt-dev libarmadillo-dev libace-dev python-h5py \
    python-matplotlib python-libxml2 gcc-multilib python-psutil \
    libgtest-dev liblapack-dev libblas-dev

#MKL
#If you want MKL instead intead of OpenBLAS, use the code below
#COPY l_mkl_11.2.1.133.tgz /opt/
#COPY mkl_silent.cfg /opt/

#The USER=root is needed for the installer to work
#ENV USER=root
#RUN cd /opt && \
#  tar -xzf l_mkl_11.2.1.133.tgz && \
#  cd l_mkl_11.2.1.133 && \
#  ./install.sh -s ../mkl_silent.cfg && \
#  cd ../ && \
#  rm -rf l_mkl*.tgz && \
#  rm -rf mkl_silent.cfg


#OpenBLAS with OpenMP
RUN cd /opt && \
    mkdir debsource && \
    cd debsource && \
    apt-get --no-install-recommends --no-install-suggests --yes build-dep libopenblas-base && \
    apt-get install --no-install-recommends --no-install-suggests --yes build-essential fakeroot devscripts && \
    apt-get source libopenblas-base && \
    cd openblas-0.2.8/ && \
    sed -i "s/NO_WARMUP=1/NO_WARMUP=1 OPENMP=1/g" debian/rules && \
    debchange -i "Compiling with OpenMP support" && \
    debuild -us -uc -i -I && \
    debi && \
    update-alternatives --set libblas.so.3 /usr/lib/openblas-base/libblas.so.3 && \
    cd /opt && \
    rm -rf debsource

#Python stuff
RUN pip install pyxb

#Set more environment variables in preparation for Gadgetron installation
ENV GADGETRON_HOME=/usr/local \
    ISMRMRD_HOME=/usr/local

ENV PATH=$PATH:$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib:/opt/intel/mkl/lib/intel64:/opt/intel/lib/intel64

#ISMRMRD
RUN mkdir /opt/code
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd.git && \
    cd ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make && \
    make install

#GADGETRON
RUN cd /opt/code && \
    git clone https://github.com/gadgetron/gadgetron.git && \
    cd gadgetron && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make && \
    make install && \
    cp ${GADGETRON_HOME}/share/gadgetron/config/gadgetron.xml.example ${GADGETRON_HOME}/share/gadgetron/config/gadgetron.xml

#ISMRMRD PYTHON API
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd-python.git &&  \
    cd ismrmrd-python && \
    sudo python setup.py install

#ISMRMRD PYTHON TOOLS
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd-python-tools.git &&  \
    cd ismrmrd-python-tools && \
    sudo python setup.py install


#SIEMENS_TO_ISMRMRD
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/siemens_to_ismrmrd.git && \
    cd siemens_to_ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make && \
    make install

#PHILIPS_TO_ISMRMRD
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/philips_to_ismrmrd.git && \
    cd philips_to_ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make && \
    make install

RUN apt-get install --yes --no-install-recommends --no-install-suggests python-twisted-web

# Clean up packages.
RUN  apt-get clean && \
   rm -rf /var/lib/apt/lists/*

#Create script for running web app
RUN cd /opt && \
    echo '#!/bin/bash' > start_web.sh && \
    echo '[ -n "${GADGETRON_RELAY_HOST}" ] && { sed -i "s/relay_host=.*/relay_host=${GADGETRON_RELAY_HOST}/g" /usr/local/share/gadgetron/config/gadgetron_web_app.cfg; }' >> start_web.sh && \
    echo 'python /usr/local/bin/gadgetron_web_app.py /usr/local/share/gadgetron/config/gadgetron_web_app.cfg' >> start_web.sh && \
    chmod +x start_web.sh

#Run the gadgetron web app
CMD /opt/start_web.sh
