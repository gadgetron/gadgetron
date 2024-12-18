# Building Pingvin

The Pingvin source code is [available on GitHub](https://github.com/gadgetron/pingvin).
The Pingvin maintainers build in a [conda](https://conda.io) environment to ensure that we are explicit about which versions of dependencies we use, etc.
You can follow the same workflow:

## Development Environment

### Using a devcontainer

The repository has a `.devcontainer` configuration for use with the [VS Code Remote](https://code.visualstudio.com/docs/remote/remote-overview). Open the repository in VS Code and make sure you have the [Remote-Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) extension installed. The extension should prompt you to reopen in a devcontainer, which will have all the dependences and tools needed for Pingvin development.

### Manual environment configuration

If you choose not to use a devcontainer, you'll need to configure your build environment manually.

The Pingvin project uses [conda](https://conda.io) to manage dependencies. To set up a build environment, make sure you have [conda installed](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) and then:

```bash
conda env create -f environment.yml
```

Then activate the environment with:

```bash
conda activate pingvin
```

and you are ready to work with the Pingvin codebase.

## Building in conda environment

In the conda environment (see above), you can build with:

```bash
git clone https://github.com/gadgetron/pingvin.git
cd pingvin
mkdir -p build
cd build
cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DUSE_MKL=ON -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX} ../
ninja
ninja install
```

## Running end-to-end tests

After building and installing Pingvin, it is recommended to run the end-to-end tests:

```bash
cd test/e2e
pytest --echo-log-on-failure
```

## Building Docker images

The Pingvin project also uses Docker images for building and distributing Pingvin. The [Dockerfile](Dockerfile) at the root of the repo provides multiple build stages for the different types (`dev` or `rt`) and flavors (`cuda` or `nocuda`) of images that you may want to use. Use the [build-images.sh](docker/build-images.sh) convenience script to build the images, e.g.:

```bash
./docker/build-images --type dev --flavor cuda
```

for a CUDA development image. See `./docker/build-images.sh --help` for details.

The `dev` (development) images contain all the dependencies needed to build Pingvin but not the actual Pingvin binaries.
The `rt` (runtime) images contain a working version of Pingvin but without the dependencies needed to build.

To run Pingvin:

```bash
docker run --gpus=all -ti ghcr.io/gadgetron/pingvin/ubuntu22.04_rt_cuda:latest --help
```

This will run the GPU enabled version of Pingvin.
