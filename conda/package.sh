#!/bin/bash
set -euo pipefail

usage()
{
  cat << EOF

Builds the pingvin conda package.

Usage: $0
EOF
}

conda_path=$(dirname "$0")
output_path="${conda_path}/build_pkg"

# Build up channel directives
channels=(
  nvidia/label/cuda-12.3.0
  ismrmrd
  conda-forge
  bioconda
  defaults
)

channel_directives=$(printf -- "-c %s " "${channels[@]}")

python3 "${conda_path}"/validate_versions.py

mkdir -p "$output_path"
bash -c "conda build --no-anaconda-upload --output-folder $output_path $channel_directives ${conda_path}"
