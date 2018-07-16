FROM ubuntu:18.04

RUN apt-get update --quiet && \
    apt-get install --no-install-recommends --no-install-suggests --yes \
    openssh-server openssh-client

EXPOSE 22

COPY entrypoint.sh /
CMD ["/entrypoint.sh"]