#! /bin/bash

set -eu

usage()
{
  cat << EOF

Builds Gadgetron runtime and dev images

Usage: $0 [options]

Options:
  --type <type>                    Type of image to build: 'dev' (development) or 'rt' (runtime) or 'all' (default)
  --flavor <flavor>                Flavor: 'cuda' or 'nocuda' or 'all' (default)
  --base-name <name>               Base name for images, <base-name>_<type>_<flavor>:<tag>, default 'ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu'
  --tag <tag>                      Tag for images, default 'latest'
  --push                           Push images
  -h, --help                       Brings up this menu
EOF
}

export DOCKER_BUILDKIT=1

types=()
types_default=("dev" "rt")

flavors=()
flavors_default=("cuda" "nocuda")

image_tag="latest"
base_name="ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu"

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
    build_stage="gadgetron_${t}_${f}"
    docker build --build-arg BUILDKIT_INLINE_CACHE=1 --target "$build_stage" -t "$image_name" -f "$(dirname "$0")/Dockerfile" "$(dirname "$0")"
    if [[ -n "${push:-}" ]]; then
      docker push "$image_name"
    fi
  done
done