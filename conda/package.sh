#!/bin/bash
set -euo pipefail

usage()
{
  cat << EOF

Builds the gadgetron conda package.

Usage: $0
EOF
}

output_path="$(dirname "$0")/build_pkg"

# Build up channel directives
channels=(
  nvidia/label/cuda-11.6.1
  ismrmrd
  gadgetron
  conda-forge
  bioconda
  defaults
  intel
)

channel_directives=$(printf -- "-c %s " "${channels[@]}")

mkdir -p "$output_path"
bash -c "conda config --set solver libmamba"
bash -c "conda build --no-anaconda-upload --output-folder $output_path $channel_directives $(dirname "$0")"
