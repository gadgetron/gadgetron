#!/bin/bash

thisDir=$(dirname $0)
IMAGE_NAME="gadgetron/codespaces-ubuntu2004"

if [ ! -f ~/.docker/config.json ]; then
    echo "Please login to docker hub with 'docker login'"
    exit 1
fi

dockerCreds=$(cat ~/.docker/config.json | jq -r '.auths."https://index.docker.io/v1/".auth')
if [ "null" == "$dockerCreds" ]; then
    echo "No credentials stored for docker hub. Please sign in with 'docker login'"
    exit 1
fi

# Make a devcontainer.json file without the image digest
cat $thisDir/devcontainer.json | jq -r ".image=\"$IMAGE_NAME\"" > $thisDir/devcontainer-strip.json

USERNAME=$(cat $thisDir/devcontainer.json | jq -r ".remoteUser")
IMAGE_ID=$(docker build --build-arg USERNAME=$USERNAME -q -t $IMAGE_NAME $thisDir)

if [ -z "$IMAGE_ID" ]; then
    echo "Image build failed"
    exit 1
fi

docker push $IMAGE_NAME
DIGEST=$(docker inspect $IMAGE_ID | jq -r .[0].RepoDigests[0])
devcontainerJson=$(cat $thisDir/devcontainer.json | jq ".image=\"${DIGEST}\"")
echo $devcontainerJson | jq . > $thisDir/devcontainer.json

# Remove temporary file
rm $thisDir/devcontainer-strip.json