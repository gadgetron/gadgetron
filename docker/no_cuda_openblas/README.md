This Docker image is a Gadgetron installation *without* CUDA support. It uses OpenBLAS (with OpenMP support) for accelerated matrix operations. To run it use:

    docker run -e "GADGETRON_RELAY_HOST=<MY RELAY HOST IP/NAME>" --name gt1 -p 9002:9002 --rm -t gadgetron

You will now have a container named `gt1` running and the Gadgetron inside will be exposed on port 9002 on the Docker host.
