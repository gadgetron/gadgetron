ARG BASE_IMAGE=ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_2004_cuda11_cudnn8_base
FROM ${BASE_IMAGE}

# Options for setup script
ARG INSTALL_ZSH="false"
ARG UPGRADE_PACKAGES="false"
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID
ARG VSCODE_DEV_CONTAINERS_SCRIPT_LIBRARY_VERSION=v0.209.3

ENV CMAKE_GENERATOR=Ninja

# Install needed packages and setup non-root user.
RUN script=$(curl -fsSL "https://raw.githubusercontent.com/microsoft/vscode-dev-containers/${VSCODE_DEV_CONTAINERS_SCRIPT_LIBRARY_VERSION}/script-library/common-debian.sh") && bash -c "$script" -- "${INSTALL_ZSH}" "${USERNAME}" "${USER_UID}" "${USER_GID}" "${UPGRADE_PACKAGES}" "true" "true"

# Download Tini for cleaning up zombie processes and default signal handling. We will copy into the runtime images
# (which do not have curl installed)
ARG TINI_VERSION=v0.19.0
RUN mkdir -p /opt/code \
    && curl -fsSL https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini --output /opt/code/tini \
    && chmod +x /opt/code/tini

ARG VSCODE_DEV_CONTAINERS_SCRIPT_LIBRARY_VERSION=v0.179.0

# Enable non-root Docker access in container
ARG ENABLE_NONROOT_DOCKER="true"
# Use the OSS Moby CLI instead of the licensed Docker CLI
ARG USE_MOBY="false"

RUN script=$(curl -fsSL "https://raw.githubusercontent.com/microsoft/vscode-dev-containers/${VSCODE_DEV_CONTAINERS_SCRIPT_LIBRARY_VERSION}/script-library/docker-debian.sh") && bash -c "$script" -- "${ENABLE_NONROOT_DOCKER}" "/var/run/docker-host.sock" "/var/run/docker.sock" "${USERNAME}" "${USE_MOBY}"

# Create a kits file for the VSCode CMake Tools extension, so you are not prompted for which kit to select whenever you open VSCode
RUN mkdir -p /home/vscode/.local/share/CMakeTools \
    && echo '[{"name":"GCC-CONDA","compilers":{"C":"/opt/conda/envs/gadgetron/bin/x86_64-conda_cos6-linux-gnu-gcc","CXX":"/opt/conda/envs/gadgetron/bin/x86_64-conda_cos6-linux-gnu-g++"}}]' > /home/vscode/.local/share/CMakeTools/cmake-tools-kits.json \
    && chown vscode /home/vscode/.local/share/CMakeTools/cmake-tools-kits.json

RUN update-locale LANG=en_US.UTF-8

# Setting the ENTRYPOINT to docker-init.sh will configure non-root access to
# the Docker socket if "overrideCommand": false is set in devcontainer.json.
# The script will also execute CMD if you need to alter startup behaviors.
ENTRYPOINT [ "/usr/local/share/docker-init.sh" ]
CMD [ "sleep", "infinity" ]
