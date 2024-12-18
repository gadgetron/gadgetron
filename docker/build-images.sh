#! /bin/bash

set -eu

REPO_ROOT="$(readlink -f "$(dirname "${BASH_SOURCE[0]}")/../")"

usage()
{
  cat << EOF

Builds Pingvin runtime and dev images

Usage: $0 [options]

Options:
  --type <type>                    Type of image to build: 'dev' (development), 'build', or 'rt' (runtime) or 'all' (default)
  --flavor <flavor>                Flavor: 'cuda', 'nocuda', or 'all' (default)
  --base-name <name>               Base name for images, <base-name>_<type>_<flavor>:<tag>, default 'ghcr.io/gadgetron/pingvin/ubuntu22.04'
  --tag <tag>                      Tag for images, default 'latest'
  --push                           Push images
  -h, --help                       Brings up this menu
EOF
}

export DOCKER_BUILDKIT=1

types=()
types_default=("dev" "build" "rt")

flavors=()
flavors_default=("cuda" "nocuda")

image_tag="latest"
base_name="ghcr.io/gadgetron/pingvin/ubuntu22.04"

while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    --flavor)
      flavors+=("$2")
      shift 2
      ;;
    --type)
      types+=("$2")
      shift 2
      ;;
    --tag)
      image_tag="$2"
      shift 2
      ;;
    --base-name)
      base_name="$2"
      shift 2
      ;;
    --push)
      push=1
      shift 1
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

if [ ${#types[@]} -eq 0 ]; then
  types=("${types_default[@]}")
fi

if [ ${#flavors[@]} -eq 0 ]; then
  flavors=("${flavors_default[@]}")
fi

for t in "${types[@]}"; do
  for f in "${flavors[@]}"; do
    image_name="${base_name}_${t}_${f}:${image_tag}"
    build_stage="pingvin_${t}_${f}"
    docker build --build-arg BUILDKIT_INLINE_CACHE=1 --cache-from type=registry,ref="$image_name" --target "$build_stage" -t "$image_name" -f "${REPO_ROOT}/Dockerfile" "${REPO_ROOT}"
    if [[ -n "${push:-}" ]]; then
      docker push "$image_name"
    fi
  done
done