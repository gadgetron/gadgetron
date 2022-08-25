#!/bin/bash

if [[ $(uname) =~ Darwin ]]; then
   set -uo pipefail  # to avoid pipeline failing on code-signing
else
   set -euo pipefail
fi

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
  nvidia
  gadgetron
  conda-forge
  bioconda
  defaults
  intel
  cefca
)

channel_directives=$(printf -- "-c %s " "${channels[@]}")

mkdir -p "$output_path"
bash -c "conda build --no-anaconda-upload --output-folder $output_path $channel_directives $(dirname "$0")"
