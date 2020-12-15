# To build image without CUDA support use base image ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_2004_base:latest
ARG BASE_IMAGE=ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_2004_cuda11_cudnn8_base:latest
FROM $BASE_IMAGE

WORKDIR /opt
RUN mkdir -p /opt/code/gadgetron
COPY . /opt/code/gadgetron/

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

CMD ["/opt/start_supervisor"]