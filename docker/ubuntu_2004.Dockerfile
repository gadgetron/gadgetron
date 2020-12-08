FROM ubuntu:20.04

WORKDIR /opt
RUN mkdir -p /opt/code/gadgetron
COPY . /opt/code/gadgetron/

# Install dependencies
RUN chmod +x /opt/code/gadgetron/docker/install_ubuntu_dependencies.sh
RUN /opt/code/gadgetron/docker/install_ubuntu_dependencies.sh

# for embedded python plot, we need agg backend
RUN mkdir -p /root/.config/matplotlib && touch /root/.config/matplotlib/matplotlibrc && echo "backend : agg" >> /root/.config/matplotlib/matplotlibrc

#Set more environment variables in preparation for Gadgetron installation
ENV GADGETRON_HOME=/usr/local \
    ISMRMRD_HOME=/usr/local

ENV PATH=$PATH:$GADGETRON_HOME/bin:$ISMRMRD_HOME/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ISMRMRD_HOME/lib:$GADGETRON_HOME/lib

RUN chmod +x /opt/code/gadgetron/docker/build_gadgetron_dependencies.sh
RUN /opt/code/gadgetron/docker/build_gadgetron_dependencies.sh

#GADGETRON
RUN cd /opt/code/gadgetron && \
    mkdir build && \
    cd build && \
    cmake ../ -G Ninja && \
    ninja && \
    ninja install && \
    cp /opt/code/gadgetron/docker/start_supervisor /opt/ && \
    cp /opt/code/gadgetron/docker/supervisord.conf /opt/

RUN pip3 install gadgetron

# Clean up packages.
RUN  apt-get clean && \
   rm -rf /var/lib/apt/lists/*

CMD ["/opt/start_supervisor"]