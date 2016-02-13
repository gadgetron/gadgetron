FROM gadgetron_ubuntu1404_cuda75_base

RUN pip install buildbot_slave
RUN mkdir /data; mkdir /data/buildbot
RUN cd /data/buildbot; buildslave create-slave slave 174.129.2.39:9989 docker-ubuntu-1404 gadgetron

RUN echo 'Michael S. Hansen (michael.hansen@nih.gov)' > /data/buildbot/slave/info/admin
RUN echo 'Ubuntu 14.04 (Docker)' > /data/buildbot/slave/info/host

#Generate supervisor file
RUN echo '[supervisord]' > /etc/supervisor/conf.d/supervisord.conf && \
    echo 'nodaemon=true' >> /etc/supervisor/conf.d/supervisord.conf && \
    echo '' >> /etc/supervisor/conf.d/supervisord.conf && \
    echo '[program:buildslave]' >> /etc/supervisor/conf.d/supervisord.conf && \
    echo 'command=twistd --nodaemon --no_save -y buildbot.tac' >> /etc/supervisor/conf.d/supervisord.conf && \
    echo 'directory=/data/buildbot/slave' >> /etc/supervisor/conf.d/supervisord.conf && \
    echo '' >> /etc/supervisor/conf.d/supervisord.conf 

CMD ["/usr/bin/supervisord"]

