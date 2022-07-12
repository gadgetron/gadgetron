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
global_channels=$(cat "$(dirname "$0")"/global.yml | yq '.channels | join(" -c ")' | tr -d '"')
channel_directive="-c $global_channels"

mkdir -p "$output_path"
bash -c "conda build --no-anaconda-upload --output-folder $output_path $channel_directive $(dirname "$0")"
