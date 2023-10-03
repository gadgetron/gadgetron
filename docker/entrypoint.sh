#!/bin/bash

set -eo pipefail

# enable conda for this shell
# shellcheck disable=SC1091
. /opt/conda/etc/profile.d/conda.sh
. /opt/set_matlab_paths.sh

# activate the environment
conda activate gadgetron

exec gadgetron "$@"
