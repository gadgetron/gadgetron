# Gadgetron Image Reconstruction Framework

The Gadgetron is an open source project for medical image reconstruction. If you find the Gadgetron useful in your research, please cite this paper:

Hansen MS, Sørensen TS. Gadgetron: An Open Source Framework for Medical Image Reconstruction. Magn Reson Med. 2013 Jun;69(6):1768-76.

Documentation for the project is available at [https://gadgetron.readthedocs.io](https://gadgetron.readthedocs.io)

## License

The Gadgetron is available under a modified MIT license. Please read [LICENSE](LICENSE) file for licensing details.


##

To build Gadgetron using the Clang compiler, from inside the 'build' directory, run:

`
cmake -GNinja -DCMAKE_C_COMPILER=${CONDA_PREFIX}/bin/clang -DCMAKE_CXX_COMPILER=${CONDA_PREFIX}/bin/clang++ -DCMAKE_BUILD_TYPE=Release -DUSE_MKL=ON -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX} ../
`

then use the ninja command to execute the build.

