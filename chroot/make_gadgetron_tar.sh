#!/bin/bash

if [ $# -eq 3 ]; then

 sudo apt-get install debootstrap -y

 sudo debootstrap --variant=buildd --arch amd64 trusty "${2}" http://gb.archive.ubuntu.com/ubuntu/

 "${3}/generate_gadgetron_root" "${1}" "${2}"

 # *** Uncomment if you have MKL compiled ***

 #sudo rm -rf "${2}/opt/intel/mkl/lib/intel64"
 #sudo rm -rf "${2}/opt/intel/lib/intel64"
 #sudo cp -rf /opt/intel/mkl/lib/intel64 "${2}/opt/intel/mkl/lib/intel64"
 #sudo cp -rf /opt/intel/lib/intel64 "${2}/opt/intel/lib/intel64"

 sudo cp "${3}/start-gadgetron.sh" "${2}"
 sudo cp "${3}/start.sh" "${2}/../"
 sudo cp "${3}/stop.sh" "${2}/../"

 sudo tar -zcvpf "${3}/gadgetron.tar.gz" --directory "${2}/../" --exclude="${2}/proc" --exclude="${2}/sys" --exclude="${2}/dev" "${2}/../"

 exit 0

else
 echo -e "\nUsage:  $0 (gadgetron_root) (chrootdir) (scripts dir)\n"
 exit 1
fi

