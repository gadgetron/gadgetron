# For base image without CUDA use BASE_IMAGE ubuntu:20.04
ARG BASE_IMAGE=nvidia/cuda:11.0-devel-ubuntu20.04
FROM $BASE_IMAGE

WORKDIR /opt
COPY install_ubuntu_dependencies.sh /opt/
RUN chmod +x /opt/install_ubuntu_dependencies.sh
COPY build_gadgetron_dependencies.sh /opt/
RUN chmod +x /opt/build_gadgetron_dependencies.sh

# Install dependencies
RUN /opt/install_ubuntu_dependencies.sh

# for embedded python plot, we need agg backend
RUN mkdir -p /root/.config/matplotlib && touch /root/.config/matplotlib/matplotlibrc && echo "backend : agg" >> /root/.config/matplotlib/matplotlibrc

#Set more environment variables in preparation for Gadgetron installation
ENV GADGETRON_HOME=/usr/local \
    ISMRMRD_HOME=/usr/local

ENV PATH=$PATH:$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib

RUN /opt/build_gadgetron_dependencies.sh


# Clean up packages.
RUN  apt-get clean && \
   rm -rf /var/lib/apt/lists/*
