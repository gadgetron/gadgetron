
To build the Docker image, move to the folder with the desired configuration and run Docker build:

    cd cuda_331_openblas/
    docker build -t gadgetron .

To start Docker container from image:

    export CUDA_DEVICES="--device=/dev/nvidia3:/dev/nvidia3 --device=/dev/nvidia2:/dev/nvidia2 --device=/dev/nvidia1:/dev/nvidia1 --device=/dev/nvidia0:/dev/nvidia0 --device=/dev/nvidiactl:/dev/nvidiactl --device=/dev/nvidia-uvm:/dev/nvidia-uvm"
    docker run -e "GADGETRON_RELAY_HOST=grenada.nhlbi.nih.gov"  ${CUDA_DEVICES} --name gt1 -p 9902:9002 --rm -t gadgetron

You will now have a container named `gt1` running and the Gadgetron inside will be exposed on port 9902 on the Docker host
