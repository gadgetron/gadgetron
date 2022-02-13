# Gadgetron Image Reconstruction Framework

## Setting up a development environment

The Gadgetron project uses [conda](ttps://conda.io) to manage dependencies. To set up a build environment, make sure you have [conda installed](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) and then:

```bash
conda env create -f environment.yml
```

Then activate the environment with:

```bash
conda activate gadgetron
```

and you are ready to work with the Gadgetron codebase.

## Building in conda environment

In the conda environment (see above), you can build with:

```bash
mkdir -p build
cd build
cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DUSE_MKL=ON -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX} ../
ninja
ninja install
```

## Running end to end tests

After building and installing the Gadgetron, it is recommended to run the end to end tests:

```bash
cd test/integration
python get_data.py
python run_tests.py --echo-log-on-failure --timeout=600 -F --stats stats.csv cases/*
```

## Building and running Docker images

The Gadgetron project uses Docker images for building and distributing the Gadgetron. The [Dockerfile](Dockerfile) at the root of the repo provides multiple build stages for the different types (`dev` or `rt`) and flavors (`cuda` or `nocuda`) of images that you may want to use. Use the [build-images.sh](build-images.sh) convenience script to build the images, e.g.:

```bash
./build-images --type dev --flavor cuda
```

for a CUDA development image. See `./build-images.sh --help` for details. 

The images provided should look like:

```
ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_dev_cuda:latest
ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_dev_nocuda:latest
ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_rt_cuda:latest
ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_rt_nocuda:latest
```

The `dev` (development) images contain all the dependencies needed to build the Gadgetron but not the actual Gadgetron binaries. The `rt` (runtime) images contain a working version of the Gadgetron but without the dependencies needed to build.

To run the Gadgetron:

```bash
docker run --gpus=all -ti -p 9004:9002 ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_rt_cuda:latest
```

This will run the GPU enabled version of the Gadgetron and expose it on port `9004`.

## Using a devcontainer for development

The repository has a `.devcontainer` configuration for use with the [VS Code Remote](https://code.visualstudio.com/docs/remote/remote-overview). Open the Gadgetron folder in VS Code and make sure you have the [Remote-Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) extension installed. The extension should prompt you to reopen in a devcontainer, which will have all the dependences and tools needed for Gadgetron development.

## Installing in conda environment

The Gadgetron can be installed in a [conda](https://conda.io) environment. To install the Gadgetron define and `environment.yaml` file with:

```yaml
name: gadgetron
channels:
  - nvidia/label/cuda-11.6.0
  - gadgetron
  - conda-forge
  - bioconda
  - defaults
  - intel
dependencies:
  - gadgetron>=4.1.2
  - siemens_to_ismrmrd>=1.0.0
```

And create the environment with:

```bash
conda env create -f environment.yaml
```

After activating the environment (with `conda activate gadgetron`), you should be able to check that everything is working with `gadgetron --info`

## Running the Gadgetron in Kubernetes

The Docker images can be deployed in a Kubernetes cluster. See [this repository](https://github.com/Microsoft/gadgetron-azure) for details.

## License

The Gadgetron is available under a modified MIT license. Please read [LICENSE](LICENSE) file for licensing details.
