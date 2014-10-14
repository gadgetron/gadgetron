Running Gadgetron in chroot
============================


Introduction
------------


Creating a chroot environment
-----------------------------

First we need to install the tools required to create chroot environments::

  sudo apt-get install dchroot debootstrap

Next we need to add an appropriate configuration to `/etc/schroot/schroot.conf'::
  
  [trusty]
  description=trusty
  location=/var/chroot/trusty
  priority=3
  users=doko
  groups=sbuild
  root-groups=root

Create the folder where we will be making the root file system::
  
  sudo mkdir -p /var/chroot/trusty


Now generate a basic root file system::

  sudo debootstrap --variant=buildd --arch amd64 trusty /var/chroot/trusty http://archive.ubuntu.com/ubuntu/


