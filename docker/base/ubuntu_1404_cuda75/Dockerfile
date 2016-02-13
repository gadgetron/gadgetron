FROM gadgetron_ubuntu1404_base

#UNPACK THE CUDA INSTALLER
ENV CUDA_DOWNLOAD_PATH http://developer.download.nvidia.com/compute/cuda/7.5/Prod/local_installers
ENV CUDA_RUN cuda_7.5.18_linux.run
ENV CUDA_INSTALL cuda-linux64-rel-7.5.18-19867135.run

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
LABEL com.nvidia.cuda.version="7.5"
RUN echo "/usr/local/cuda/lib" >> /etc/ld.so.conf.d/cuda.conf && \
    echo "/usr/local/cuda/lib64" >> /etc/ld.so.conf.d/cuda.conf && \
    ldconfig

RUN echo "/usr/local/nvidia/lib" >> /etc/ld.so.conf.d/nvidia.conf && \
    echo "/usr/local/nvidia/lib64" >> /etc/ld.so.conf.d/nvidia.conf

ENV PATH /usr/local/cuda/bin:${PATH}
ENV LD_LIBRARY_PATH /usr/local/nvidia/lib:/usr/local/nvidia/lib64:${LD_LIBRARY_PATH}

# Clean up packages.
RUN  apt-get clean && \
   rm -rf /var/lib/apt/lists/*
