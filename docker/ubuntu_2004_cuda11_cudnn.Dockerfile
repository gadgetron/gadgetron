FROM nvidia/cuda:11.0-devel-ubuntu20.04

WORKDIR /opt

# Install dependencies
COPY install_ubuntu_dependencies.sh ./
RUN chmod +x install_ubuntu_dependencies.sh
RUN ./install_ubuntu_dependencies.sh

# for embedded python plot, we need agg backend
RUN mkdir -p /root/.config/matplotlib && touch /root/.config/matplotlib/matplotlibrc && echo "backend : agg" >> /root/.config/matplotlib/matplotlibrc

#Set more environment variables in preparation for Gadgetron installation
ENV GADGETRON_HOME=/usr/local \
    ISMRMRD_HOME=/usr/local

ENV PATH=$PATH:$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib

ARG GADGETRON_URL=https://github.com/gadgetron/gadgetron
ARG GADGETRON_BRANCH=master

COPY build_gadgetron.sh ./
RUN chmod +x build_gadgetron.sh
RUN ./build_gadgetron.sh --source-repo-url ${GADGETRON_URL} --source-branch ${GADGETRON_BRANCH}

# Clean up packages.
RUN  apt-get clean && \
   rm -rf /var/lib/apt/lists/*

CMD ["/opt/start_supervisor"]