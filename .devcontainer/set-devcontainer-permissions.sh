#!/bin/bash

CONDA_FOLDER="/opt/conda/envs/$(grep 'name:' /tmp/build/filtered_environment.yml | awk '{print $2}')"

# Make sure any user can follow links:
sudo find "${CONDA_FOLDER}/pkgs" -type d ! -perm o=x -exec chmod o+x {} \;

# Make sure we can install binaries (share/, lib/, include/, bin) 
sudo find "${CONDA_FOLDER}/include" -type d ! -perm o=w -exec chmod o+w {} \;
sudo find "${CONDA_FOLDER}/lib" -type d ! -perm o=w -exec chmod o+w {} \;
sudo find "${CONDA_FOLDER}/bin" -type d ! -perm o=w -exec chmod o+w {} \;
sudo find "${CONDA_FOLDER}/share" -type d ! -perm o=w -exec chmod o+w {} \;

