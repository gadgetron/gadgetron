This Docker image is a Gadgetron installation *without* CUDA support. It uses OpenBLAS (with OpenMP support) for accelerated matrix operations. To run it use:

    docker run -e "GADGETRON_RELAY_HOST=<MY RELAY HOST IP/NAME>" --name gt1 -p 9002:9002 --rm -t gadgetron

Other possible ports to expose are

* 9080: ReST API
* 22: SSH server
* 8002: CloudBus relay

You will now have a container named `gt1` running and the Gadgetron inside will be exposed on port 9002 on the Docker host.
