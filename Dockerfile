# Shared arguments
ARG USERNAME="vscode"
ARG USER_UID=1000
ARG USER_GID=$USER_UID

FROM ghcr.io/gadgetron/gadgetron/ubuntu22.04 AS pingvin_baseimage
LABEL org.opencontainers.image.source=https://github.com/gadgetron/pingvin

ARG USERNAME
ARG USER_UID
ARG USER_GID
ARG HOME=/home/$USERNAME

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y sudo wget git-core rsync curl net-tools libxml2 \
    && apt-get clean

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME -s /bin/bash \
    #
    # [Optional] Add sudo support. Omit if you don't need to install software after connecting.
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

ARG MINIFORGE_VERSION=24.3.0-0
ARG CONDA_GID=900

# Based on https://github.com/conda-forge/miniforge-images/blob/master/ubuntu/Dockerfile
RUN wget --no-hsts --quiet https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-${MINIFORGE_VERSION}-Linux-$(uname -m).sh -O /tmp/miniforge.sh \
    && /bin/bash /tmp/miniforge.sh -b -p /opt/conda \
    && rm /tmp/miniforge.sh \
    && /opt/conda/bin/conda clean --tarballs --index-cache --packages --yes \
    && find /opt/conda -follow -type f -name '*.a' -delete \
    && find /opt/conda -follow -type f -name '*.pyc' -delete \
    && /opt/conda/bin/conda clean --force-pkgs-dirs --all --yes  \
    && groupadd -r conda --gid ${CONDA_GID} \
    && usermod -aG conda ${USERNAME} \
    && chown -R :conda /opt/conda \
    && chmod -R g+w /opt/conda \
    && find /opt -type d | xargs -n 1 chmod g+s

# Copy environment, which will be filtered for later staged
COPY --chown=$USER_UID:conda environment.yml /tmp/build/

# Create mount points for tests
RUN mkdir -p /test && chown ${USER_UID}:${USER_GID} /test
VOLUME /test

# Copy bashrc used to activate conda environment
COPY --chown=$USER_UID:conda .devcontainer/devcontainer.bashrc /opt/devcontainer/

# Add a section to /etc/bash.bashrc that ensures that a section is present at the end of ~/.bashrc.
# We can't just write to .bashrc from here because it will be overwritten if the vscode user has
# opted to use their own dotfiles repo. The dotfiles repo is cloned after the postCreateCommand
# in the devcontainer.json file is executed.
RUN echo "\n\
if ! grep -q \"^source /opt/devcontainer/devcontainer.bashrc\" \${HOME}/.bashrc; then\n\
	echo \"source /opt/devcontainer/devcontainer.bashrc\" >> \${HOME}/.bashrc\n\
fi\n" >> /etc/bash.bashrc

ENV TINI_VERSION=v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini

# Install the yardl tool
ARG YARDL_VERSION=0.6.2
RUN wget --quiet "https://github.com/microsoft/yardl/releases/download/v${YARDL_VERSION}/yardl_${YARDL_VERSION}_linux_x86_64.tar.gz" \
    && tar -xzf "yardl_${YARDL_VERSION}_linux_x86_64.tar.gz" \
    && mv yardl "/usr/local/bin/" \
    && rm "yardl_${YARDL_VERSION}_linux_x86_64.tar.gz"

FROM pingvin_baseimage AS pingvin_dev_cuda
ARG USER_UID
ARG HOME
USER ${USER_UID}
RUN mkdir -p ${HOME}/.cache/conda/notices && sudo chown -R ${USER_UID}:conda ${HOME}/.cache/conda/notices
RUN grep -v "#.*\<NOFILTER\>" /tmp/build/environment.yml > /tmp/build/filtered_environment.yml
RUN umask 0002 && /opt/conda/bin/conda env create -f /tmp/build/filtered_environment.yml && /opt/conda/bin/conda clean -afy && sudo chown -R :conda /opt/conda
USER root

FROM pingvin_baseimage AS pingvin_dev_nocuda
ARG USER_UID
ARG HOME
USER ${USER_UID}
RUN mkdir -p ${HOME}/.cache/conda/notices && sudo chown -R ${USER_UID}:conda ${HOME}/.cache/conda/notices
RUN grep -v "#.*\<cuda\>" /tmp/build/environment.yml > /tmp/build/filtered_environment.yml
RUN umask 0002 && /opt/conda/bin/conda env create -f /tmp/build/filtered_environment.yml && /opt/conda/bin/conda clean -afy && sudo chown -R :conda /opt/conda
USER root

FROM pingvin_dev_cuda AS pingvin_build_cuda
ARG USER_UID
USER ${USER_UID}
WORKDIR /opt
RUN sudo chown $USER_UID:$USER_GID /opt && mkdir -p /opt/code/pingvin && mkdir -p /opt/package
COPY --chown=$USER_UID:conda . /opt/code/pingvin/
SHELL ["/bin/bash", "-c"]
RUN . /opt/conda/etc/profile.d/conda.sh && umask 0002 && conda activate pingvin && sh -x && \
    cd /opt/code/pingvin && \
    mkdir build && \
    cd build && \
    cmake ../ -GNinja -DUSE_MKL=ON -DCMAKE_INSTALL_PREFIX=/opt/package && \
    ninja && \
    ninja install

FROM pingvin_dev_nocuda AS pingvin_build_nocuda
ARG USER_UID
USER ${USER_UID}
WORKDIR /opt
RUN sudo chown $USER_UID:$USER_GID /opt && mkdir -p /opt/code/pingvin && mkdir -p /opt/package
COPY --chown=$USER_UID:conda . /opt/code/pingvin/
SHELL ["/bin/bash", "-c"]
RUN . /opt/conda/etc/profile.d/conda.sh && umask 0002 && conda activate pingvin && sh -x && \
    cd /opt/code/pingvin && \
    mkdir build && \
    cd build && \
    cmake ../ -GNinja -DUSE_MKL=ON -DCMAKE_INSTALL_PREFIX=/opt/package && \
    ninja && \
    ninja install

FROM pingvin_baseimage AS pingvin_rt_cuda
ARG USER_UID
ARG HOME
USER ${USER_UID}
RUN mkdir -p ${HOME}/.cache/conda/notices && sudo chown -R ${USER_UID}:conda ${HOME}/.cache/conda/notices
RUN grep -v "#.*\<dev\>" /tmp/build/environment.yml > /tmp/build/filtered_environment.yml
RUN umask 0002 && /opt/conda/bin/conda env create -f /tmp/build/filtered_environment.yml && /opt/conda/bin/conda clean -afy && sudo chown -R :conda /opt/conda
COPY --from=pingvin_build_cuda --chown=$USER_UID:conda /opt/package /opt/conda/envs/pingvin/
COPY --from=pingvin_build_cuda --chown=$USER_UID:conda /opt/code/pingvin/docker/entrypoint.sh /opt/
RUN chmod +x /opt/entrypoint.sh
RUN sudo mkdir -p /opt/e2e-test && sudo chown ${USER_GID}:${USER_UID} /opt/e2e-test
COPY --from=pingvin_build_cuda --chown=$USER_UID:conda /opt/code/pingvin/test/e2e /opt/e2e-test/
ENTRYPOINT [ "/tini", "--", "/opt/entrypoint.sh" ]

FROM pingvin_baseimage AS pingvin_rt_nocuda
ARG USER_UID
ARG HOME
USER ${USER_UID}
RUN mkdir -p ${HOME}/.cache/conda/notices && sudo chown -R ${USER_UID}:conda ${HOME}/.cache/conda/notices
RUN grep -v "#.*\<cuda\|dev\>" /tmp/build/environment.yml > /tmp/build/filtered_environment.yml
RUN umask 0002 && /opt/conda/bin/conda env create -f /tmp/build/filtered_environment.yml && /opt/conda/bin/conda clean -afy && sudo chown -R :conda /opt/conda
COPY --from=pingvin_build_nocuda --chown=$USER_UID:conda /opt/package /opt/conda/envs/pingvin/
COPY --from=pingvin_build_nocuda --chown=$USER_UID:conda /opt/code/pingvin/docker/entrypoint.sh /opt/
RUN chmod +x /opt/entrypoint.sh
RUN sudo mkdir -p /opt/e2e-test && sudo chown ${USER_GID}:${USER_UID} /opt/e2e-test
COPY --from=pingvin_build_nocuda --chown=$USER_UID:conda /opt/code/pingvin/test/e2e /opt/e2e-test/
ENTRYPOINT [ "/tini", "--", "/opt/entrypoint.sh" ]
