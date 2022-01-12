#!/bin/bash

set -eu

REPO_ROOT="$(readlink -f "$(dirname "${BASH_SOURCE[0]}")/../")"

usage()
{
  cat << EOF

Test a Gadgetron docker container image

Usage: $0 [options]

Options:
  --image <image name>              Name of image to test
  --cases <cases glob>              Test cases
  --ignore-requirements | -i <name> Ignore requirements parameters for test, e.g. python,cuda
  --gpus                            gpus argument for docker run
  -h, --help                        Brings up this menu
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
    --ignore-requirements | -i)
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

if [[ -z "${ignore_requirements:-}" ]]; then
  ignore=""
else
  ignore="--ignore-requirements=${ignore_requirements}"
fi

# Since the Gadgetron instance is not running as root, we need some gymnastics to get the test stats out of the container. 
volume_name="gadgetron_test_${RANDOM}"
docker volume create "$volume_name" 1>/dev/null
docker run --rm --gpus="${gpus}" -v ${HOST_WORKSPACE_DIR}/test/integration/data:/opt/integration-test/data \
    --mount src="$volume_name",destination=/test "$image_name" \
    /bin/bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate gadgetron && cd /opt/integration-test/ && python run_tests.py --echo-log-on-failure ${ignore} --timeout=600 -F --stats /test/stats.csv ${cases} && python stats_to_junit.py -i /test/stats.csv -o /test/junit.xml"
docker run --rm --mount src="$volume_name",destination=/test busybox cat /test/stats.csv > stats.csv
docker run --rm --mount src="$volume_name",destination=/test busybox cat /test/junit.xml > junit.xml
docker volume rm "$volume_name" 1>/dev/null
