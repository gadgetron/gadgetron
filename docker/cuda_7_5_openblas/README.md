This Docker image is a full Gadgetron installation with CUDA support (driver version 331 and CUDA 5.5) and using OpenBLAS (with OpenMP support) for accelerated matrix operations. To run it use:

    export CUDA_DEVICES="--device=/dev/nvidia0:/dev/nvidia0 --device=/dev/nvidiactl:/dev/nvidiactl --device=/dev/nvidia-uvm:/dev/nvidia-uvm"
    docker run -e "GADGETRON_RELAY_HOST=<MY RELAY HOST IP/NAME>"  ${CUDA_DEVICES} --name gt1 -p 9002:9002 --rm -t gadgetron

You will now have a container named `gt1` running and the Gadgetron inside will be exposed on port 9002 on the Docker host. The `CUDA_DEVICES` part is only needed if you would like to expose your CUDA devices inside the Docker container. There is a helper script `nvidia_devices.sh`, which will generate this string for you. 
