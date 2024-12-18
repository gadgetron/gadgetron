#!/bin/bash

set -eo pipefail

# enable conda for this shell
# shellcheck disable=SC1091
. /opt/conda/etc/profile.d/conda.sh

# activate the environment
conda activate pingvin

exec pingvin "$@"