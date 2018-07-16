#!/bin/bash

#Set new keys
rm -rf /etc/ssh/ssh_host_*
dpkg-reconfigure openssh-server

#Switch off DNS checking
sed -i.bak "s/#UseDNS no/UseDNS no/g" /etc/ssh/sshd_config

#Import public key
mkdir -p /root/.ssh
touch /root/.ssh/authorized_keys
echo ${PUBLIC_KEY} > /root/.ssh/authorized_keys
chmod 600 /root/.ssh/authorized_keys

#Run SSH server
mkdir -p /run/sshd
/usr/sbin/sshd -D