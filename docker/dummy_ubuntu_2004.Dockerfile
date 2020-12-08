FROM ubuntu:20.04

RUN mkdir -p /opt/code/gadgetron/build/test && \
    echo '#!/bin/bash' >> /opt/code/gadgetron/build/test/test_all && \
    echo "echo 'Dummy tests'" >> /opt/code/gadgetron/build/test/test_all && \
    echo "exit 0" >> /opt/code/gadgetron/build/test/test_all && \
    chmod +x /opt/code/gadgetron/build/test/test_all


