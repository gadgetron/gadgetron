FROM ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_2004_cuda11_cudnn8_base

# Options for setup script
ARG INSTALL_ZSH="false"
ARG UPGRADE_PACKAGES="false"
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID
ARG VSCODE_DEV_CONTAINERS_SCRIPT_LIBRARY_VERSION=v0.209.3

# Install needed packages and setup non-root user.
RUN script=$(curl -fsSL "https://raw.githubusercontent.com/microsoft/vscode-dev-containers/${VSCODE_DEV_CONTAINERS_SCRIPT_LIBRARY_VERSION}/script-library/common-debian.sh") && bash -c "$script" -- "${INSTALL_ZSH}" "${USERNAME}" "${USER_UID}" "${USER_GID}" "${UPGRADE_PACKAGES}" "true" "true" \
    && apt-get install -y gdb

RUN update-locale LANG=en_US.UTF-8