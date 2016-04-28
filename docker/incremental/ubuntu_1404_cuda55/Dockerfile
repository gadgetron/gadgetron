FROM gadgetron/ubuntu1404_cuda55_base

#ISMRMRD
RUN mkdir /opt/code
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd.git && \
    cd ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j $(nproc) && \
    make install

#GADGETRON
RUN cd /opt/code && \
    git clone https://github.com/gadgetron/gadgetron.git && \
    cd gadgetron && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j $(nproc) && \
    make install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.gadgetron.sha1 --value `git rev-parse HEAD` && \
    cp ${GADGETRON_HOME}/share/gadgetron/config/gadgetron.xml.example ${GADGETRON_HOME}/share/gadgetron/config/gadgetron.xml

#HASH for ISMRMRD
RUN cd /opt/code/ismrmrd && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.ismrmrd.sha1 --value `git rev-parse HEAD` 

#ISMRMRD PYTHON API
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd-python.git &&  \
    cd ismrmrd-python && \
    sudo python setup.py install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.ismrmrdpython.sha1 --value `git rev-parse HEAD` 


#ISMRMRD PYTHON TOOLS
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd-python-tools.git &&  \
    cd ismrmrd-python-tools && \
    sudo python setup.py install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.ismrmrdpythontools.sha1 --value `git rev-parse HEAD` 

#SIEMENS_TO_ISMRMRD
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/siemens_to_ismrmrd.git && \
    cd siemens_to_ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j $(nproc) && \
    make install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.siemens_to_ismrmrd.sha1 --value `git rev-parse HEAD` 


#PHILIPS_TO_ISMRMRD
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/philips_to_ismrmrd.git && \
    cd philips_to_ismrmrd && \
    mkdir build && \
    cd build && \
    cmake ../ && \
    make -j $(nproc) && \
    make install && \
    /opt/code/gadgetron/docker/manifest --key .io.gadgetron.philips_to_ismrmrd.sha1 --value `git rev-parse HEAD` 


#Create script for running web app
RUN cd /opt && \
    echo '#!/bin/bash' > start_web.sh && \
    echo '[ -n "${GADGETRON_RELAY_HOST}" ] && { sed -i "s/relay_host=.*/relay_host=${GADGETRON_RELAY_HOST}/g" /usr/local/share/gadgetron/config/gadgetron_web_app.cfg; }' >> start_web.sh && \
    echo '[ -n "${GADGETRON_LB_ENDPOINT}" ] && { sed -i "s/lb_endpoint=.*/lb_endpoint=${GADGETRON_LB_ENDPOINT}/g" /usr/local/share/gadgetron/config/gadgetron_web_app.cfg; }' >> start_web.sh && \
    echo '[ -n "${GADGETRON_LOG_FILE}" ] && { sed -i "s|logfile=/tmp/gadgetron.log|logfile=${GADGETRON_LOG_FILE}|g" /usr/local/share/gadgetron/config/gadgetron_web_app.cfg; }' >> start_web.sh && \
    echo 'python /usr/local/bin/gadgetron_web_app.py /usr/local/share/gadgetron/config/gadgetron_web_app.cfg' >> start_web.sh && \
    chmod +x start_web.sh

#Generate supervisor file
RUN echo '[supervisord]' > /etc/supervisor/conf.d/supervisord.conf && \
    echo 'nodaemon=true' >> /etc/supervisor/conf.d/supervisord.conf && \
    echo '' >> /etc/supervisor/conf.d/supervisord.conf && \
    echo '[program:gadgetron]' >> /etc/supervisor/conf.d/supervisord.conf && \
    echo 'command=/opt/start_web.sh' >> /etc/supervisor/conf.d/supervisord.conf && \
    echo '' >> /etc/supervisor/conf.d/supervisord.conf 

# Clean up packages.
RUN  apt-get clean && \
   rm -rf /var/lib/apt/lists/*

CMD ["/usr/bin/supervisord"]

