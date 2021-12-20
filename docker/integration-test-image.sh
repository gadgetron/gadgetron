#!/bin/bash

set -eu

IMAGE_NAME="$1"
REPO_ROOT="$(readlink -f "$(dirname "${BASH_SOURCE[0]}")")"

usage()
{
  cat << EOF

Test a Gadgetron docker container image

Usage: $0 [options]

Options:
  --image <image name>             Type of image to build: 'dev' (development) or 'rt' (runtime) or 'all' (default)
  --cases <cases glob>             Flavor: 'cuda' or 'nocuda' or 'all' (default)
  --ignore-requirement | -i <name> ignore requirements parameters for test
  --gpus                           gpus argument for docker run
  -h, --help                       Brings up this menu
EOF
}

image_name="ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_rt_cuda:latest"
gpus=all
cases="cases/*"

while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    --image)
      image_name="$2"
      shift 2
      ;;
    --cases)
      cases="$2"
      shift 2
      ;;
    --ignore-requirement | -i)
      ignore_requirements="$2"
      shift 2
      ;;
    --gpus)
      gpus="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit
      ;;
    *)
      echo "ERROR: unknown option \"$key\""
      usage
      exit 1
      ;;
  esac
done

# If we are not in a devcontainer, assume repo root is the root
if [[ -z "${HOST_WORKSPACE_DIR:-}" ]]; then
    HOST_WORKSPACE_DIR="$REPO_ROOT"
fi

# Since the Gadgetron instance is not running as root, we need some gymnastics to get the test stats out of the container. 
volume_name="gadgetron_test_${RANDOM}"
docker volume create "$volume_name"
docker run --rm -it --gpus="${gpus}" -v ${HOST_WORKSPACE_DIR}/test/integration/data:/opt/integration-test/data \
    --mount src="$volume_name",destination=/opt/integration-test/results "$image_name" \
    /bin/bash -c "ls -al /opt/integration-test/"
# docker run --rm -it --gpus="${gpus}" -v ${HOST_WORKSPACE_DIR}/test/integration/data:/opt/integration-test/data \
#     --mount src="$volume_name",destination=/opt/integration-test/results "$image_name" \
#     /bin/bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate gadgetron && cd /opt/integration-test/ && python run_tests.py --echo-log-on-failure --timeout=600 -F --stats /opt/integration-test/results/stats.csv ${cases}"
TMP_CID=$(docker run --rm -d --mount src="$volume_name",destination=/test busybox true)
docker cp "$TMP_CID":/test/stats.csv .
docker rm "$TMP_CID"
docker volume rm "$volume_name"