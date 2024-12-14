# Using Pingvin

In the following walkthrough we demonstrate how to send raw data through Pingvin and obtain images.


The walkthrough assumes that you have installed or built Pingvin. You can also use the Pingvin docker images
that are distributed from the GitHub container registries. Specifically, a command like:

```bash
docker run -it --rm --gpus=all --entrypoint /bin/bash ghcr.io/gadgetron/pingvin/ubuntu22.04_rt_cuda:latest
```

should open a terminal where you can type the commands in this walkthrough.

## Validate installation

First, validate that Pingvin is installed and working. The command `pingvin --info` should give you information
about your installed version of Pingvin and it would look something like this:

```
$ pingvin --info
Pingvin Version Info
  -- Version            : 5.0.0.alpha.1
  -- Git SHA1           : NA
  -- System Memory size : 112705 MB
  -- CUDA Support       : YES
  -- NVCC Flags         : -gencode arch=compute_60,code=sm_60;-gencode arch=compute_61,code=sm_61;-gencode arch=compute_70,code=sm_70;-gencode arch=compute_75,code=sm_75;-gencode arch=compute_80,code=sm_80;-gencode arch=compute_86,code=sm_86 --std=c++17
    * Number of CUDA capable devices: 1
      - Device 0: Tesla P100-PCIE-16GB
         + CUDA Driver Version / Runtime Version: 11.6/11.6
         + CUDA Capability Major / Minor version number: 6.0
         + Total amount of global GPU memory: 16280 MB
```

The output may vary on your specific setup, but you will see error messages if Pingvin is not installed or not installed correctly.


## Performing Reconstructions

To perform a reconstruction using Pingvin, you must provide serialized data either from an
input file or through an [Anonymous pipe](https://en.wikipedia.org/wiki/Anonymous_pipe).


### Generating test data

The MRD library has an application that allows you to generate a simple dataset for testing purposes. If you have
your own dataset this step is not neccesary.

To generate test data type the following command in a terminal:

```bash
mrd_phantom -r 10 > testdata.mrd
```

This generates a dataset in the current working directory (./testdata.mrd) with 8 coils and 10 repetitions. You can
then use this data to perform reconstructions with Pingvin in either server or stream mode.


### Running Pingvin

To run Pingvin, you must specify the configuration. We can use the default configuration to reconstruct the dataset we just generated:

```bash
pingvin -c default.xml -i testdata.mrd -o reconstructed.mrd
```

You should see Pingvin log its output:

```bash
$ pingvin -c default.xml -i testdata.mrd -o reconstructed.mrd
12-02 15:31:12.079 DEBUG [system_info.cpp:225] Executable path: "/opt/conda/envs/pingvin/bin/pingvin"
12-02 15:31:12.079 DEBUG [system_info.cpp:231] Default Pingvin home: "/opt/conda/envs/pingvin"
12-02 15:31:12.079 WARNING [initialization.cpp:38] Environment variable 'OMP_WAIT_POLICY' not set to 'PASSIVE'.
12-02 15:31:12.079 WARNING [initialization.cpp:39] Gadgetron may experience serious performance issues under heavy load (multiple simultaneous reconstructions, etc.)
12-02 15:31:12.079 INFO [main.cpp:85] Pingvin 5.0.0.alpha.1 [c0baaa8d864fcba8efd562b7b35fd963956a78b5]
12-02 15:31:12.079 INFO [StreamConsumer.h:68] Loading configuration from: /opt/conda/envs/pingvin/share/pingvin/config/default.xml
```

You can convert the reconstructed dataset if you'd like to use HDF5 tools for viewing the data:

```bash
cat reconstructed.mrd | mrd_stream_to_hdf5 > reconstructed.h5
```

### More Streaming Examples

#### Docker image

The Pingvin docker image can be leveraged to invoke a reconstruction without building/installing Pingvin locally.

```bash
mrd_phantom -r 10 | docker run -i --gpus=all ghcr.io/gadgetron/pingvin/ubuntu22.04_rt_cuda:latest -c default.xml -o reconstructed.mrd
```

#### Chainging Pingvin

A complex reconstruction can be formulated by chaining compatible reconstruction configurations, for
instance a chain which produces complex images can be converted to float images, and then to short images later:

```bash
mrd_phantom | pingvin -c Generic_Cartesian_Grappa_Complex.xml | pingvin -c stream_complex_to_float.xml | pingvin -c stream_float_to_short.xml -o reconstructed.mrd
```
