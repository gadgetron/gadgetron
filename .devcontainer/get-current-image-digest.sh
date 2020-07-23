#!/bin/bash

docker inspect $(docker images -q cloudenvimage:latest) | jq -r .[0].RepoDigests[0]
