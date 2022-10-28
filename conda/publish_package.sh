#!/bin/bash
set -euo pipefail

usage()
{
  cat << EOF

Publishes a conda package.

Usage: $0 [options]

Options:
  -p|--package_path <path>               Path to the package (tar.gz) to push
  -u|--user <user>                       Anaconda.org channeluser or organization
  -t|--token <token>                     Token for uploading to anaconda.org
  -f|--force                             Force push even if package exists
  -h| --help  Brings up this menu
EOF
}

while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -p|--package_path)
      package_path="$2"
      shift
      shift
      ;;
    -u|--user)
      user="$2"
      shift
      shift
      ;;
    -t|--token)
      token="$2"
      shift
      shift
      ;;
    --force)
      force=1
      shift
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

if [[ -z "${package_path:-}" ]]; then
  echo "You cannot push to anaconda without a package"
  echo "Please supply a package path with the --package_path argument"
  exit 1
fi
if [[ -z "${token:-}" ]]; then
  echo "You cannot push to anaconda without a token"
  echo "Please supply a token with the --token argument"
  exit 1
fi
if [[ -z "${user:-}" ]]; then
  echo "You cannot push to anaconda without a user"
  echo "Please supply a user with the --user argument"
  exit 1
fi

force_directive="--skip-existing"
if [[ -n ${force:-} ]]; then
  force_directive="--force"
fi

anaconda -t "$token" upload -u "$user" $force_directive "$package_path"
