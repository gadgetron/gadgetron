#!/bin/bash

set -eo pipefail

# enable conda for this shell
# shellcheck disable=SC1091
. /opt/conda/etc/profile.d/conda.sh

# activate the environment
conda activate gadgetron

# exec the cmd/command in this process, making it pid 1
exec gadgetron "$@"