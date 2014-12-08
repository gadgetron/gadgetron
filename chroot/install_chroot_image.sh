#!/bin/bash

if [ $(id -u) -ne 0 ]; then 
 	echo -e "\nPlease start the script as a root or sudo!\n"
 	exit 1

else

	BASEDIR=$(dirname $0)

	install_img=0

	if [ $# -eq 1 ]; then
		CHROOT_IMAGE_FILENAME=${1}
		CHROOT_INSTALL_PATH=/home/gadgetron_chroot
	else
		if [ $# -eq 2 ]; then
			CHROOT_IMAGE_FILENAME=${1}
			CHROOT_INSTALL_PATH=${2}
		else
			if [ $# -eq 3 ]; then
				if [ ${2} == "latest" ]; then
					TAR_NAME=`find ${1} -type f -name 'gadgetron-*.tar.gz' |sort |head -n1`
					CHROOT_IMAGE_FILENAME=${TAR_NAME}
				else
					CHROOT_IMAGE_FILENAME=${1}/${2}			
				fi
				CHROOT_INSTALL_PATH=${3}
			else
				if [ $# -eq 4 ]; then
					if [ ${2} == "latest" ]; then
						TAR_NAME=`find ${1} -type f -name 'gadgetron-*.tar.gz' |sort |head -n1`
						CHROOT_IMAGE_FILENAME=${1}/${TAR_NAME}
					else
						CHROOT_IMAGE_FILENAME=${1}/${2}			
					fi
					CHROOT_INSTALL_PATH=${3}

					if [ ${4} -eq 1 ]; then
						install_img=1
						IMG_NAME=`find ${1} -type f -name 'gadgetron-*.img' |sort |head -n1`
						CHROOT_IMAGE_IMG_FILENAME=${IMG_NAME}
					fi
				else
					echo -e "\nUsage 1, install chroot image to /home/gadgetron_chroot: $0 chroot_image_file chroot_install_path"
				  	echo -e "\nUsage 2, install chroot image to selected install path: $0 chroot_image_file chroot_install_path"
				  	echo -e "\nUsage 3, : $0 chroot_image_path chroot_image_name chroot_install_path"
				  	echo -e "\n           install chroot image to selected install path, if chroot_image_name=latest, the newest chroot image in the folder will be installed: $0 chroot_image_file chroot_install_path"
				  	echo -e "\nUsage 4, : $0 chroot_image_path chroot_image_name chroot_install_path install_img"
				  	echo -e "\n           like Usage 3, if install_img=1, the corresponding .img package will be copied to chroot_install_path"
				  	exit 1
				fi  
			fi  
		fi  
	fi

  	service gadgetron_chroot stop

	echo CHROOT_IMAGE_FILENAME=${CHROOT_IMAGE_FILENAME}
	echo CHROOT_INSTALL_PATH=${CHROOT_INSTALL_PATH}

	mkdir -p ${CHROOT_INSTALL_PATH}

	cp -rf ${CHROOT_IMAGE_FILENAME} ${CHROOT_INSTALL_PATH}/

	FILENAME_WITH_EXTENSION=${CHROOT_IMAGE_FILENAME##*/}
	FILENAME=${FILENAME_WITH_EXTENSION%.*}
	FILENAME=${FILENAME%.*}
	echo ${FILENAME}

	mkdir ${CHROOT_INSTALL_PATH}/${FILENAME}

	echo untar ${CHROOT_INSTALL_PATH}/${FILENAME_WITH_EXTENSION} ... 

	tar -xzf ${CHROOT_INSTALL_PATH}/${FILENAME_WITH_EXTENSION} --directory="${CHROOT_INSTALL_PATH}/${FILENAME}" .

	rm -f ${CHROOT_INSTALL_PATH}/current

	ln -s ${CHROOT_INSTALL_PATH}/${FILENAME} ${CHROOT_INSTALL_PATH}/current

	cp -f ${CHROOT_INSTALL_PATH}/current/chroot-root/gadgetron/webapp/gadgetron_chroot.conf /etc/init/

	if [ ${install_img} -eq 1 ]; then
                echo "copy image file : ${CHROOT_IMAGE_IMG_FILENAME} ... "		
		cp -f ${CHROOT_IMAGE_IMG_FILENAME} ${CHROOT_INSTALL_PATH}/
	fi

	service gadgetron_chroot start

	exit 0
fi
