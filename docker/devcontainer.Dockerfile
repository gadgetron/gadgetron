FROM mcr.microsoft.com/vscode/devcontainers/base:0.201.8-focal AS devcontainer

# Install needed packages and setup non-root user.
ARG USERNAME="vscode"
ARG USER_UID=1000
ARG USER_GID=$USER_UID
ARG HOME=/home/$USERNAME
ARG CONDA_ENVIRONMENT_NAME=gadgetron

# The version of conda to use
ARG CONDA_VERSION=4.10.3

# Based on https://github.com/ContinuumIO/docker-images/blob/master/miniconda3/debian/Dockerfile.
# We also install conda-lock.
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh \
    && mkdir -p /opt \
    && sh miniconda.sh -b -p /opt/conda \
    && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && [ -z "$CONDA_VERSION" ] || /opt/conda/bin/conda install -n base conda=$CONDA_VERSION \
    && /opt/conda/bin/conda install -c conda-forge -n base conda-lock \
    && /opt/conda/bin/conda clean -afy \
	&& chown -R $USER_UID:$USER_GID /opt/conda

# Add a section to /etc/bash.bashrc that ensures that a section is present at the end of ~/.bashrc.
# We can't just write to .bashrc from here because it will be overwritten if the vscode user has
# opted to use their own dotfiles repo. The dotfiles repo is cloned after the postCreateCommand
# in the devcontainer.json file is executed.
RUN echo "\n\
if ! grep -q \"^source /opt/conda/etc/profile.d/conda.sh\" ${HOME}/.bashrc; then\n\
	echo \"source /opt/conda/etc/profile.d/conda.sh\" >> ${HOME}/.bashrc\n\
	echo \"conda activate ${CONDA_ENVIRONMENT_NAME}\" >> ${HOME}/.bashrc\n\
    echo \"source <(kubectl completion bash)\" >> ${HOME}/.bashrc\n\
fi\n" >> /etc/bash.bashrc

# Create a conda environment from the lockfile in the repo root.
COPY environment.yml /tmp/build/
RUN /opt/conda/bin/conda env create -f /tmp/build/environment.yml \
    && chown -R $USER_UID:$USER_GID /opt/conda/envs \
    && chown -R $USER_UID:$USER_GID ${HOME}/.conda

ENV CMAKE_GENERATOR=Ninja

# Install dependencies not evailable through conda or pip
COPY docker/bootstrap-conda.sh /tmp/build/
RUN chmod +x /tmp/build/bootstrap-conda.sh
ENV PATH="/app:/opt/conda/condabin:${PATH}"
RUN conda run --no-capture-output -n ${CONDA_ENVIRONMENT_NAME} /tmp/build/bootstrap-conda.sh

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

# Setting the ENTRYPOINT to docker-init.sh will configure non-root access to
# the Docker socket if "overrideCommand": false is set in devcontainer.json.
# The script will also execute CMD if you need to alter startup behaviors.
ENTRYPOINT [ "/usr/local/share/docker-init.sh" ]
CMD [ "sleep", "infinity" ]

# Create a kits file for the VSCode CMake Tools extension, so you are not prompted for which kit to select whenever you open VSCode
# RUN mkdir -p /home/vscode/.local/share/CMakeTools \
#    && echo '[{"name":"GCC-10","compilers":{"C":"/opt/conda/envs/wabakimi/bin/x86_64-conda_cos6-linux-gnu-gcc","CXX":"/opt/conda/envs/wabakimi/bin/x86_64-conda_cos6-linux-gnu-g++"}}]' > /home/vscode/.local/share/CMakeTools/cmake-tools-kits.json \
#    && chown vscode /home/vscode/.local/share/CMakeTools/cmake-tools-kits.json
