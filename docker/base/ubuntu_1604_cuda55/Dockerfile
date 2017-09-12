FROM gadgetron/ubuntu1604_base

LABEL com.nvidia.volumes.needed="nvidia_driver"

#UNPACK THE CUDA INSTALLER
ENV CUDA_DOWNLOAD_PATH http://developer.download.nvidia.com/compute/cuda/5_5/rel/installers 
ENV CUDA_RUN cuda_5.5.22_linux_64.run
ENV CUDA_INSTALL cuda-linux64-rel-5.5.22-16488124.run
ENV CUDA_VERSION 5.5

RUN cd /opt && \
  wget ${CUDA_DOWNLOAD_PATH}/${CUDA_RUN} --no-verbose && \
  chmod +x *.run && \
  mkdir nvidia_installers && \
  ./${CUDA_RUN} -silent -extract=`pwd`/nvidia_installers && \
  rm ${CUDA_RUN}

#INSTALL CUDA LIBRARIES, etc. 
RUN cd /opt/nvidia_installers && \
    ./${CUDA_INSTALL} -silent -noprompt && \
    cd ../ && \
    rm -rf nvidia_installers

#SET SOME ENV variables for CUDA
LABEL com.nvidia.cuda.version="5.5"
RUN echo "/usr/local/cuda/lib" >> /etc/ld.so.conf.d/cuda.conf && \
    echo "/usr/local/cuda/lib64" >> /etc/ld.so.conf.d/cuda.conf && \
    ldconfig

RUN echo "/usr/local/nvidia/lib" >> /etc/ld.so.conf.d/nvidia.conf && \
    echo "/usr/local/nvidia/lib64" >> /etc/ld.so.conf.d/nvidia.conf

ENV PATH /usr/local/nvidia/bin:/usr/local/cuda/bin:${PATH}
ENV LD_LIBRARY_PATH /usr/local/nvidia/lib:/usr/local/nvidia/lib64:${LD_LIBRARY_PATH}

# Clean up packages.
RUN  apt-get clean && \
   rm -rf /var/lib/apt/lists/*
