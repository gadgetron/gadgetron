This Docker image is a full Gadgetron installation with CUDA support (CUDA 7.5) and using OpenBLAS (with OpenMP support) for accelerated matrix operations. To start Docker container from image use the `nvidia-docker` script, which can be downloaded from https//raw.githubusercontent.com/NVIDIA/nvidia-docker/master/nvidia-docker or https//raw.githubusercontent.com/gadgetron/gadgetron/master/docker/nvidia-docker 

    GPU=0 nvidia-docker run -e "GADGETRON_RELAY_HOST=<MY RELAY HOST IP/NAME>" --name gt1 -p 9002:9002 --rm -t gadgetron

Other possible ports to expose are

* 9080: ReST API
* 22: SSH server
* 8002: CloudBus relay
 
You will now have a container named `gt1` running and the Gadgetron inside will be exposed on port 9002 on the Docker host. The use of the nvidia-docker script is only needed if you would like to expose GPUs inside the container and use CUDA with the Gadgetron. 
