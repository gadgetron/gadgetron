GADGETRON IMAGE RECONSTRUCTION FRAMEWORK
=========================================

Detailed installation instructions and manual is available at:

https://github.com/gadgetron/gadgetron/wiki

General Building Instructions (on Unix platforms) to install Gadgetron to system path
-------------------------------------------------

```
mkdir build
cd build
cmake ../
make
sudo make install
```

Please check manual for detailed instructions for your platform.

To install Gadgetron in your home folder
-------------------------------------------------
Please follow the instructions at [Set up Gadgetron in your home folder](https://github.com/gadgetron/gadgetron/wiki/Visual-debug-Gadgetron-in-Ubuntu-using-Eclipse)

After compiling Gadgetron, it is a good idea to run integration test
-------------------------------------------------
Please follow instructionst at [Integration test in Gadgetron](https://github.com/gadgetron/gadgetron/wiki/Integration-test-in-Gadgetron)

To deploy Gadgetron, it is convenient to use docker
-------------------------------------------------
Please follow instructionst at [Deploy Gadgetron via Docker](docker/)

Gadgetron does support cluster/cloud deployment
-------------------------------------------------
If you have multiple computers, they can be assembled as a cluster to run distributed Gadgetron. Please check instructions [Set up Gadgetron Cluster in your LAN](https://github.com/gadgetron/gadgetron/wiki/How-to-set-up-Gadgetron-cloud-in-the-LAN-environment)

Please read LICENSE file for licensing details.
