#!/bin/bash

if [ $# -eq 3 ]; then

 sudo apt-get install debootstrap -y

 sudo mkdir -p "${2}/gadgetron"

 sudo debootstrap --variant=buildd --arch amd64 trusty "${2}/gadgetron" http://gb.archive.ubuntu.com/ubuntu/

 "${3}/generate_gadgetron_root" "${1}" "${2}/gadgetron"

 # *** Uncomment if you have MKL compiled ***

 #sudo rm -rf "${2}/gadgetron/opt/intel/mkl/lib/intel64"
 #sudo rm -rf "${2}/gadgetron/opt/intel/lib/intel64"
 #sudo cp -rf /opt/intel/mkl/lib/intel64 "${2}/gadgetron/opt/intel/mkl/lib/intel64"
 #sudo cp -rf /opt/intel/lib/intel64 "${2}/gadgetron/opt/intel/lib/intel64"

 sudo cp "${3}/start-gadgetron.sh" "${2}/gadgetron"
 sudo cp "${3}/start.sh" "${2}"
 sudo cp "${3}/stop.sh" "${2}"

 sudo tar -zcvf "${3}/gadgetron_chroot.tar.gz" --directory "${2}" .
 exit 0

else
 echo -e "\nUsage:      $0 (gadgetron_root) (chrootdir) (scripts dir)\nImportant:  Make sure that the chrootdir is empty!\n"
 exit 1
fi

