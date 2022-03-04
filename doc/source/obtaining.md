# Obtaining and running the Gadgetron

We recommend two different ways of obtaining and running the Gadgetron; using a conda environment or using Docker containers.

## Installing in conda environment

The Gadgetron can be installed in a [conda](https://conda.io) environment. To install the Gadgetron define and `environment.yml` file with:

```yaml
name: gadgetron
channels:
  - nvidia/label/cuda-11.6.1
  - gadgetron
  - conda-forge
  - bioconda
  - defaults
  - intel
dependencies:
  - gadgetron>=4.1.5
  - siemens_to_ismrmrd>=1.2
```

And create the environment with:

```bash
conda env create -f environment.yml
```

After activating the environment (with `conda activate gadgetron`), you should be able to check that everything is working with `gadgetron --info`

## Using Docker containers

Docker images are built automatically from the Gadgetron project. The latest runtime images are:

1. `ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_rt_cuda:latest`: The latest runtime image with CUDA support.
1. `ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_rt_nocuda:latest`: The latest runtime image without CUDA support. 

To run the Gadgetron:

```bash
docker run --gpus=all -ti -p 9004:9002 ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_rt_cuda:latest
```

This will run the GPU enabled version of the Gadgetron and expose it on port `9004`.

For details on how to build these images yourself, see out [build instructions](building)

## Running the Gadgetron in Kubernetes

The Docker images can be deployed in a Kubernetes cluster. See [this repository](https://github.com/Microsoft/gadgetron-azure) for details.
