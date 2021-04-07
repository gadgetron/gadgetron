#!/bin/bash

thisFolder=$(dirname $0)

if [ -d "${thisFolder}/../.vscode" ]; then
    echo ".vscode folder already exists. Will not overwrite."
else
    cp -r "${thisFolder}/.vscode" "${thisFolder}/../"
fi