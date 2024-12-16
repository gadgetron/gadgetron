# Obtaining and running Pingvin

We recommend two different ways of obtaining and running Pingvin; using a conda environment or using Docker containers.

## Installing in conda environment

Pingvin can be installed in a [conda](https://conda.io) environment. To install Pingvin, define an `environment.yml` file with:

```yaml
name: pingvin
channels:
  - ismrmrd
  - gadgetron
  - nvidia/label/cuda-12.3.0
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - pingvin>=5.0.0.alpha.1
  - siemens_to_ismrmrd>=1.2.13
```

And create the environment with:

```bash
conda env create -f environment.yml
```

After activating the environment (with `conda activate pingvin`), you should be able to check that everything is working with `pingvin --info`

## Using Docker containers

Docker images are built automatically from the Pingvin project. The latest runtime images are:

1. `ghcr.io/gadgetron/pingvin/ubuntu22.04_rt_cuda:latest`: The latest runtime image with CUDA support.
1. `ghcr.io/gadgetron/pingvin/ubuntu22.04_rt_nocuda:latest`: The latest runtime image without CUDA support.

To run Pingvin:

```bash
docker run --gpus=all -it ghcr.io/gadgetron/pingvin/ubuntu22.04_rt_cuda:latest --info
```

This will run the GPU enabled version of Pingvin.

For details on how to build these images yourself, see out [build instructions](building)

## Running Pingvin in Kubernetes

The Docker images can be deployed in a Kubernetes cluster. See [this repository](https://github.com/Microsoft/gadgetron-azure) for details.
