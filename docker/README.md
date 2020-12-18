# Gadgetron Docker Images

To build a Docker image, move to the `docker/` folder and build the desired image:

```bash
cd docker/
docker build -t gttestbuild -f ubuntu_2004.Dockerfile .
```

or

```bash
cd docker/
docker build -t gttestbuild -f ubuntu_2004_cuda11_cudnn.Dockerfile .
```


To build an image from a your own git repo and branch:

```bash
docker build --build-arg GADGETRON_URL=https://github.com/hansenms/gadgetron --build-arg GADGETRON_BRANCH=mybranch -t gttestbuild -f ubuntu_2004_cuda11_cudnn.Dockerfile -t mygadgetron --no-cache .
```

To run the Gadgetron in a container:

```bash
docker run -ti --rm --gpus all gttestbuild
```

The `--gpus all` argument is only needed if you are running with GPUs

You can also use one of the official docker images:

```bash
docker run -ti --rm --gpus all ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_2004_cuda11_cudnn8
```


 

