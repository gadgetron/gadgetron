#!/bin/bash

thisFolder=$(dirname $0)

if [ -d "${thisFolder}/../.vscode" ]; then
    echo ".vscode folder already exists. Will not overwrite."
else
    cp -r "${thisFolder}/.vscode" "${thisFolder}/../"
fi

mkdir -p /home/vscode/.local/share/CMakeTools && \
    echo '[{"name":"GCC-CONDA","compilers":{"C":"/opt/conda/envs/pingvin/bin/x86_64-conda_cos6-linux-gnu-gcc","CXX":"/opt/conda/envs/pingvin/bin/x86_64-conda_cos6-linux-gnu-g++"}}]' > /home/vscode/.local/share/CMakeTools/cmake-tools-kits.json && \
    chown vscode:conda /home/vscode/.local/share/CMakeTools/cmake-tools-kits.json