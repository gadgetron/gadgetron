# Using the Gadgetron

In the following walkthrough we demonstrate how to send raw data through the Gadgetron and obtain images.


The walkthrough assumes that you have installed or built the Gadgetron. You can also use the Gadgetron docker images
that are distributed from the GitHub container registries. Specifically, a command like:

```bash
docker run --net=host -it --rm --gpus=all ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_rt_cuda:latest /bin/bash
```

should open a terminal where you can type the commands in this walkthrough. You may require more than one terminal to
have simultaneous client and server sessions while using the Gadgetron in server mode.


## Validate installation

First, validate that the Gadgetron is installed and working. The command `gadgetron --info` should give you information
about your installed version of the Gadgetron and it would look something like this:

```
$ gadgetron --info
Gadgetron Version Info
  -- Version            : 4.1.5
  -- Git SHA1           : NA
  -- System Memory size : 112705 MB
  -- Python Support     : YES
  -- Julia Support      : NO
  -- Matlab Support     : NO
  -- CUDA Support       : YES
  -- NVCC Flags         : -gencode arch=compute_60,code=sm_60;-gencode arch=compute_61,code=sm_61;-gencode arch=compute_70,code=sm_70;-gencode arch=compute_75,code=sm_75;-gencode arch=compute_80,code=sm_80;-gencode arch=compute_86,code=sm_86 --std=c++17
    * Number of CUDA capable devices: 1
      - Device 0: Tesla P100-PCIE-16GB
         + CUDA Driver Version / Runtime Version: 11.6/11.6
         + CUDA Capability Major / Minor version number: 6.0
         + Total amount of global GPU memory: 16280 MB
```

The output may vary on your specific setup, but you will see error messages if the Gadgetron is not installed or not installed correctly.


## Performing Reconstructions

The Gadgetron can be run in two ways to perform reconstructions. One is a server mode in which the Gadgetron listens for
active connections on a desired TCP port and performs reconstructions with data received over the network.
The second is a stream mode in which the Gadgetron is provided serialized data from the local machine, either from an
input file or through an [Anonymous pipe](https://en.wikipedia.org/wiki/Anonymous_pipe).


### Generating test data

The ISMRMRD library has an application that allows you to generate a simple dataset for testing purposes. If you have
your own dataset this step is not neccesary. However note that the input file used in the following examples will
reference the data as 'testdata.h5'.

To generate test data type the following command in a terminal:

```bash
ismrmrd_generate_cartesian_shepp_logan -r 10
```

This generates a dataset in the current working directory (./testdata.h5) with 8 coils and 10 repetitions. You can
then use this data to perform reconstructions with the gadgetron in either server or stream mode.


### Server Mode

Open up two terminals. We will use one for the Gadgetron server and one for the Gadgetron client. In the server
terminal you should be able to simply type `gadgetron` and see the following output:

```
$ gadgetron
03-08 00:33:12.267 DEBUG [gadgetron_paths.cpp:107] Executable path: "/opt/conda/envs/gadgetron/bin/gadgetron"
03-08 00:33:12.267 DEBUG [gadgetron_paths.cpp:113] Default Gadgetron home: "/opt/conda/envs/gadgetron"
03-08 00:33:12.267 WARNING [initialization.cpp:38] Environment variable 'OMP_WAIT_POLICY' not set to 'PASSIVE'.
03-08 00:33:12.267 WARNING [initialization.cpp:39] Gadgetron may experience serious performance issues under heavy load (multiple simultaneous reconstructions, etc.)
03-08 00:33:12.267 INFO [main.cpp:76] Gadgetron 4.1.5 [NA]
03-08 00:33:12.267 INFO [main.cpp:77] Running on port 9002
03-08 00:33:12.731 INFO [storage.cpp:32] Running storage server on port 9112
03-08 00:33:12.732 INFO [Server.cpp:25] Gadgetron home directory: "/opt/conda/envs/gadgetron"
03-08 00:33:12.732 INFO [Server.cpp:26] Gadgetron working directory: "/tmp/gadgetron/"
Configuring services, Running on port 9002
```

This means that your Gadgetron server is running and ready to receive data.

The next step is to utilize the gadgetron_ismrmrd_client to send the data to the gadgetron over the network:

```
$ gadgetron_ismrmrd_client -f testdata.h5
Gadgetron ISMRMRD client
  -- host            :      localhost
  -- port            :      9002
  -- hdf5 file  in   :      testdata.h5
  -- hdf5 group in   :      /dataset
  -- conf            :      default.xml
  -- loop            :      1
  -- hdf5 file out   :      out.h5
  -- hdf5 group out  :      2022-03-08 00:35:04
```

This should begin a reconstruction and you should see output like this in the server terminal:

```
03-08 00:35:04.947 INFO [Server.cpp:38] Accepted connection from: ::ffff:127.0.0.1
03-08 00:35:04.948 INFO [ConfigConnection.cpp:131] Connection state: [CONFIG]
03-08 00:35:04.948 DEBUG [ConfigConnection.cpp:73] Reading config file: "/opt/conda/envs/gadgetron/share/gadgetron/config/default.xml"
03-08 00:35:04.949 INFO [HeaderConnection.cpp:84] Connection state: [HEADER]
03-08 00:35:04.949 DEBUG [RESTStorageClient.cpp:301] Using storage address: http://localhost:9112
03-08 00:35:04.949 INFO [StreamConnection.cpp:75] Connection state: [STREAM]
03-08 00:35:04.965 DEBUG [Stream.cpp:58] Loading Gadget AccTrig of class AcquisitionAccumulateTriggerGadget from gadgetron_mricore
03-08 00:35:04.965 DEBUG [Stream.cpp:58] Loading Gadget Buff of class BucketToBufferGadget from gadgetron_mricore
03-08 00:35:04.965 DEBUG [Stream.cpp:58] Loading Gadget ImageFinish of class ImageFinishGadget from gadgetron_mricore
03-08 00:35:04.965 DEBUG [Stream.cpp:58] Loading Gadget SimpleRecon of class SimpleReconGadget from gadgetron_mricore
03-08 00:35:04.966 DEBUG [Stream.cpp:58] Loading Gadget RemoveROOversampling of class RemoveROOversamplingGadget from gadgetron_mricore
03-08 00:35:04.967 DEBUG [Stream.cpp:58] Loading Gadget ImageArraySplit of class ImageArraySplitGadget from gadgetron_mricore
03-08 00:35:04.969 DEBUG [RemoveROOversamplingGadget.cpp:43] RemoveROOversamplingGadget:omp_set_num_threads(1) ...
03-08 00:35:04.969 DEBUG [Stream.cpp:58] Loading Gadget Extract of class ExtractGadget from gadgetron_mricore
03-08 00:35:05.246 DEBUG [AcquisitionAccumulateTriggerGadget.cpp:129] Trigger (1) occurred, sending out 1 buckets
03-08 00:35:05.246 DEBUG [BucketToBufferGadget.cpp:50] BUCKET_SIZE 256 ESPACE 0
03-08 00:35:05.249 DEBUG [BucketToBufferGadget.cpp:90] End of bucket reached, sending out 1 ReconData buffers
03-08 00:35:05.264 DEBUG [AcquisitionAccumulateTriggerGadget.cpp:129] Trigger (2) occurred, sending out 1 buckets
03-08 00:35:05.265 DEBUG [BucketToBufferGadget.cpp:50] BUCKET_SIZE 256 ESPACE 0
03-08 00:35:05.268 DEBUG [BucketToBufferGadget.cpp:90] End of bucket reached, sending out 1 ReconData buffers
03-08 00:35:05.338 DEBUG [AcquisitionAccumulateTriggerGadget.cpp:129] Trigger (3) occurred, sending out 1 buckets
03-08 00:35:05.338 DEBUG [BucketToBufferGadget.cpp:50] BUCKET_SIZE 256 ESPACE 0
03-08 00:35:05.341 DEBUG [BucketToBufferGadget.cpp:90] End of bucket reached, sending out 1 ReconData buffers
03-08 00:35:05.391 DEBUG [AcquisitionAccumulateTriggerGadget.cpp:129] Trigger (4) occurred, sending out 1 buckets
03-08 00:35:05.392 DEBUG [BucketToBufferGadget.cpp:50] BUCKET_SIZE 256 ESPACE 0
03-08 00:35:05.395 DEBUG [BucketToBufferGadget.cpp:90] End of bucket reached, sending out 1 ReconData buffers
03-08 00:35:05.446 DEBUG [AcquisitionAccumulateTriggerGadget.cpp:129] Trigger (5) occurred, sending out 1 buckets
03-08 00:35:05.447 DEBUG [BucketToBufferGadget.cpp:50] BUCKET_SIZE 256 ESPACE 0
03-08 00:35:05.449 DEBUG [BucketToBufferGadget.cpp:90] End of bucket reached, sending out 1 ReconData buffers
03-08 00:35:05.494 DEBUG [AcquisitionAccumulateTriggerGadget.cpp:129] Trigger (6) occurred, sending out 1 buckets
03-08 00:35:05.494 DEBUG [BucketToBufferGadget.cpp:50] BUCKET_SIZE 256 ESPACE 0
03-08 00:35:05.497 DEBUG [BucketToBufferGadget.cpp:90] End of bucket reached, sending out 1 ReconData buffers
03-08 00:35:05.567 DEBUG [AcquisitionAccumulateTriggerGadget.cpp:129] Trigger (7) occurred, sending out 1 buckets
03-08 00:35:05.567 DEBUG [BucketToBufferGadget.cpp:50] BUCKET_SIZE 256 ESPACE 0
03-08 00:35:05.570 DEBUG [BucketToBufferGadget.cpp:90] End of bucket reached, sending out 1 ReconData buffers
03-08 00:35:05.659 DEBUG [AcquisitionAccumulateTriggerGadget.cpp:129] Trigger (8) occurred, sending out 1 buckets
03-08 00:35:05.659 DEBUG [BucketToBufferGadget.cpp:50] BUCKET_SIZE 256 ESPACE 0
03-08 00:35:05.662 DEBUG [BucketToBufferGadget.cpp:90] End of bucket reached, sending out 1 ReconData buffers
03-08 00:35:05.751 DEBUG [AcquisitionAccumulateTriggerGadget.cpp:129] Trigger (9) occurred, sending out 1 buckets
03-08 00:35:05.751 DEBUG [BucketToBufferGadget.cpp:50] BUCKET_SIZE 256 ESPACE 0
03-08 00:35:05.754 DEBUG [BucketToBufferGadget.cpp:90] End of bucket reached, sending out 1 ReconData buffers
03-08 00:35:05.836 DEBUG [Gadget.h:130] Shutting down Gadget ()
03-08 00:35:05.836 DEBUG [AcquisitionAccumulateTriggerGadget.cpp:129] Trigger (10) occurred, sending out 1 buckets
03-08 00:35:05.836 DEBUG [BucketToBufferGadget.cpp:50] BUCKET_SIZE 256 ESPACE 0
03-08 00:35:05.839 DEBUG [BucketToBufferGadget.cpp:90] End of bucket reached, sending out 1 ReconData buffers
03-08 00:35:05.851 DEBUG [Gadget.h:130] Shutting down Gadget ()
03-08 00:35:05.851 DEBUG [Gadget.h:130] Shutting down Gadget ()
03-08 00:35:05.852 INFO [Core.cpp:77] Connection state: [FINISHED]
```

By default the output images are saved in the folder in which you started the `gadgetron_ismrmrd_client`.
The client appends the result to an HDF5 file called `out.h5` (if no other file name is specified).
A group is created with the current time and date and the images are stored in that group. If you run multiple
reconstructions one after another, the results will be added to the same file with a new group is created for each run.
This makes it easy to compare results from different reconstructions.

### Stream Mode

You can use the ISMRMRD HDF5/Stream conversion tools (or adapters).

Open up a single terminal. We begin by converting 'testdata.h5' generated above to an intermediary streaming format:

    ismrmrd_hdf5_to_stream -i testdata.h5 -o test_out.dat

This produces a file which can be consumed by the gadgetron in stream mode as follows:

    gadgetron --from_stream -c default.xml -i test_out.dat -o img_out.dat

The output of this command generates a file that contains images in an intermediary format, which must be converted back:

    ismrmrd_stream_to_hdf5 -i img_out.dat -o out.h5

More simply, this entire chain can be performed in a single shell command by using anonymous pipes:

    ismrmrd_stream_to_hdf5 -i testdata.h5 --use-stdout | gadgetron --from_stream -c default.xml | ismrmrd_stream_to_hdf5 --use-stdin -o out.h5

Also, a more complicated reconstruction can be formulated by chaining compatible reconstruction configurations, for
instance a chain which produces complex images can be converted to floats images, and then to short images later:

    ismrmrd_hdf5_to_stream -i testdata.h5 --use-stdout | gadgetron --from_stream -c Generic_Cartesian_Grappa_Complex.xml | gadgetron --from_stream -c stream_complex_to_float.xml | gadgetron --from_stream -c stream_float_to_short.xml | ismrmrd_stream_to_hdf5 --use-stdin -o out.h5

A chain like this can be broken up with the intermediary data being processed in different ways enabling post processing:

    ismrmrd_hdf5_to_stream --use-stdout -i testdata.h5 | gadgetron --from_stream -c Generic_Cartesian_Grappa_Complex.xml -o tmp.dat

And consumed in different ways:

    gadgetron --from_stream -c stream_complex_to_float.xml -i tmp.dat | gadgetron --from_stream -c stream_float_to_short.xml | ismrmrd_stream_to_hdf5 --use-stdin -o out1.h5
    cat tmp.dat | gadgetron --from_stream -c stream_complex_to_float.xml | gadgetron --from_stream -c stream_float_to_short.xml | ismrmrd_stream_to_hdf5 --use-stdin -o out2.h5

Additionally, the gadgetron docker container can be leveraged to invoke a reconstruction without requiring gadgetron to
be installed locally:

    ismrmrd_hdf5_to_stream --use-stdout -i testdata.h5 | docker run -i --gpus=all ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_rt_cuda:latest --from_stream -c default.xml | ismrmrd_stream_to_hdf5 --use-stdin -o out.h5

### Viewing output

If you have followed either the server or stream based reconstruction above you should have an output file out.h5.
The images can be viewed in a Jupyter Notebook using the following code snippet:

```python
import h5py
from matplotlib import pyplot as plt
import numpy as np
%matplotlib inline

f = h5py.File('/workspaces/gadgetron/out.h5')
image_series = list(f.keys())
data = np.array(f[image_series[0]]['image_0']['data']).squeeze()
plt.imshow(data[0,:,:])
```

In this case, you should see an image of the Shepp-Logan phantom.

## Storage Server

Gadgetron uses a [storage server](https://github.com/ismrmrd/mrd-storage-server) to persist data that subsequent reconstructions may depend on. The storage server address is specified as a command-line argument to gadgetron:

```bash
gadgetron --storage-address https://my-storage-server:3333
```

For convenience during development, Gadgetron will start a storage server if no address is specified. However, for
production scenarios, we stongly recommend setting up the storage server ahead of time and passing in its URL to
Gadgetron. This is especially true while running in stream mode, as it is likely that there is a desire for the storage
to persist beyond the lifetime of a single instance of the Gadgetron consuming streamed data.
