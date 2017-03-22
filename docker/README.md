
To build the Docker image, move to the folder with the desired configuration and run Docker build:

    cd incremental/ubuntu_1604_cuda75
    docker build --no-cache -t gadgetron .

To build an image from a your own git repo and branch:

   docker build --build-arg GADGETRON_URL=https://github.com/hansenms/gadgetron --build-arg GADGETRON_BRANCH=mybranch --no-cache -t mygadgetron .

The `--no-cache` option will ensure that you do a fresh pull from the git repos. To start Docker container from image use the `nvidia-docker` script, which can be downloaded from https//raw.githubusercontent.com/NVIDIA/nvidia-docker/master/nvidia-docker or https//raw.githubusercontent.com/gadgetron/gadgetron/master/docker/nvidia-docker 

    GPU=0 nvidia-docker run -e "GADGETRON_RELAY_HOST=<MY RELAY HOST IP/NAME>" --name=gt1 --publish=9002:9002 --rm -t gadgetron

You will now have a container named `gt1` running and the Gadgetron inside will be exposed on port 9002 on the Docker host. The use of the nvidia-docker script is only needed if you would like to expose GPUs inside the container and use CUDA with the Gadgetron. 

Other possible ports to expose are

* 9001: The [supervisord](https://supervisord.org) web interface.
* 9080: ReST API
* 8002: CloudBus relay

To convert a docker image to a `chroot` root file system:

     docker export $(docker create gadgetron/ubuntu_1604_no_cuda) > gadgetron_fs.tar

This tar file can then be unpacked:

     mkdir rootfs
     sudo tar -xf gadgetron_fs.tar -C rootfs/

And to run the Gadgetron in the chroot filesystem:

    sudo ${GADGETRON_SOURCE}/docker/start_gadgetron_chroot rootfs/


 

