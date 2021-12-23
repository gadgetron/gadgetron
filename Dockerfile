# Shared arguments
ARG USERNAME="vscode"
ARG USER_UID=1000
ARG USER_GID=$USER_UID

FROM ubuntu:20.04 AS gadgetron_baseimage

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

# The version of conda to use
ARG CONDA_VERSION=4.11.0

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
	&& chown -R $USER_UID:$USER_GID /opt/conda \
    && mkdir -p ${HOME}/.conda \
    && chown -R $USER_UID:$USER_GID ${HOME}/.conda

# Copy environment, which will be filtered for later staged
COPY --chown=$USER_UID:${USER_GID} environment.yml /tmp/build/

# Create mount points for tests
RUN mkdir -p /test && chown ${USER_UID}:${USER_GID} /test
VOLUME /test

# Add a section to /etc/bash.bashrc that ensures that a section is present at the end of ~/.bashrc.
# We can't just write to .bashrc from here because it will be overwritten if the vscode user has
# opted to use their own dotfiles repo. The dotfiles repo is cloned after the postCreateCommand
# in the devcontainer.json file is executed.
RUN echo "\n\
if ! grep -q \"^source /opt/conda/etc/profile.d/conda.sh\" ${HOME}/.bashrc; then\n\
	echo \"source /opt/conda/etc/profile.d/conda.sh\" >> ${HOME}/.bashrc\n\
	echo \"conda activate $(grep 'name:' /tmp/build/environment.yml | awk '{print $2}')\" >> ${HOME}/.bashrc\n\
fi\n" >> /etc/bash.bashrc

FROM gadgetron_baseimage AS gadgetron_cudadevimage_base
ARG USER_UID
ARG USER_GID
USER ${USER_UID}:${USER_GID}
RUN grep -v "#.*\<NOFILTER\>" /tmp/build/environment.yml > /tmp/build/filtered_environment.yml
RUN /opt/conda/bin/conda env create -f /tmp/build/filtered_environment.yml && /opt/conda/bin/conda clean -afy 

FROM gadgetron_cudadevimage_base AS gadgetron_dependency_build
ARG USER_UID
ARG USER_GID
USER ${USER_UID}:${USER_GID}
COPY --chown=$USER_UID:${USER_GID} docker/bootstrap-conda.sh /tmp/build/
RUN chmod +x /tmp/build/bootstrap-conda.sh
ENV PATH="/app:/opt/conda/condabin:${PATH}"
RUN conda run --no-capture-output -n "$(grep 'name:' /tmp/build/environment.yml | awk '{print $2}')" /tmp/build/bootstrap-conda.sh

FROM gadgetron_baseimage AS gadgetron_dev_nocuda
ARG USER_UID
ARG USER_GID
USER ${USER_UID}:${USER_GID}
RUN grep -v "#.*\<cuda\>" /tmp/build/environment.yml > /tmp/build/filtered_environment.yml
RUN /opt/conda/bin/conda env create -f /tmp/build/filtered_environment.yml && /opt/conda/bin/conda clean -afy
COPY --from=gadgetron_dependency_build --chown=$USER_UID:${USER_GID} /tmp/dep-build/package/ /opt/conda/envs/gadgetron/

FROM gadgetron_cudadevimage_base AS gadgetron_dev_cuda
ARG USER_UID
ARG USER_GID
COPY --from=gadgetron_dependency_build --chown=$USER_UID:${USER_GID} /tmp/dep-build/package/ /opt/conda/envs/gadgetron/

FROM gadgetron_dev_cuda AS gadgetron_cudabuild
ARG USER_UID
ARG USER_GID
USER ${USER_UID}:${USER_GID}
WORKDIR /opt
RUN sudo chown $USER_UID:$USER_GID /opt && mkdir -p /opt/code/gadgetron && mkdir -p /opt/package
COPY --chown=$USER_UID:${USER_GID} . /opt/code/gadgetron/
SHELL ["/bin/bash", "-c"]
RUN . /opt/conda/etc/profile.d/conda.sh && conda activate gadgetron && sh -x && \
    cd /opt/code/gadgetron && \
    mkdir build && \
    cd build && \
    cmake ../ -GNinja -DUSE_MKL=ON -DCMAKE_INSTALL_PREFIX=/opt/package && \
    ninja && \
    ninja install

FROM gadgetron_dev_nocuda AS gadgetron_nocudabuild
ARG USER_UID
ARG USER_GID
USER ${USER_UID}:${USER_GID}
WORKDIR /opt
RUN sudo chown $USER_UID:$USER_GID /opt && mkdir -p /opt/code/gadgetron && mkdir -p /opt/package
COPY --chown=$USER_UID:${USER_GID} . /opt/code/gadgetron/
SHELL ["/bin/bash", "-c"]
RUN . /opt/conda/etc/profile.d/conda.sh && conda activate gadgetron && sh -x && \
    cd /opt/code/gadgetron && \
    mkdir build && \
    cd build && \
    cmake ../ -GNinja -DUSE_MKL=ON -DCMAKE_INSTALL_PREFIX=/opt/package && \
    ninja && \
    ninja install

FROM gadgetron_baseimage AS gadgetron_rt_cuda
ARG USER_UID
ARG USER_GID
USER ${USER_UID}:${USER_GID}
RUN grep -v "#.*\<dev\>" /tmp/build/environment.yml > /tmp/build/filtered_environment.yml
RUN /opt/conda/bin/conda env create -f /tmp/build/filtered_environment.yml && /opt/conda/bin/conda clean -afy
COPY --from=gadgetron_dependency_build --chown=$USER_UID:${USER_GID} /tmp/dep-build/package/ /opt/conda/envs/gadgetron/
COPY --from=gadgetron_cudabuild --chown=$USER_UID:${USER_GID} /opt/package /opt/conda/envs/gadgetron/
COPY --from=gadgetron_cudabuild --chown=$USER_UID:${USER_GID} /opt/code/gadgetron/docker/start_supervisor /opt/
COPY --from=gadgetron_cudabuild --chown=$USER_UID:${USER_GID} /opt/code/gadgetron/docker/supervisord.conf /opt/
RUN sudo mkdir -p /opt/integration-test && sudo chown ${USER_GID}:${USER_UID} /opt/integration-test
COPY --from=gadgetron_cudabuild --chown=$USER_UID:${USER_GID} /opt/code/gadgetron/test/integration /opt/integration-test/
CMD ["/opt/conda/bin/conda", "run", "-n", "gadgetron", "--no-capture-output", "/opt/start_supervisor"]

FROM gadgetron_baseimage AS gadgetron_rt_nocuda
ARG USER_UID
ARG USER_GID
USER ${USER_UID}:${USER_GID}
RUN grep -v "#.*\<cuda\|dev\>" /tmp/build/environment.yml > /tmp/build/filtered_environment.yml
RUN /opt/conda/bin/conda env create -f /tmp/build/filtered_environment.yml && /opt/conda/bin/conda clean -afy
COPY --from=gadgetron_dependency_build --chown=$USER_UID:${USER_GID} /tmp/dep-build/package/ /opt/conda/envs/gadgetron/
COPY --from=gadgetron_nocudabuild --chown=$USER_UID:${USER_GID} /opt/package /opt/conda/envs/gadgetron/
COPY --from=gadgetron_nocudabuild --chown=$USER_UID:${USER_GID} /opt/code/gadgetron/docker/start_supervisor /opt/
COPY --from=gadgetron_nocudabuild --chown=$USER_UID:${USER_GID} /opt/code/gadgetron/docker/supervisord.conf /opt/
RUN sudo mkdir -p /opt/integration-test && sudo chown ${USER_GID}:${USER_UID} /opt/integration-test
COPY --from=gadgetron_nocudabuild --chown=$USER_UID:${USER_GID} /opt/code/gadgetron/test/integration /opt/integration-test/
CMD ["/opt/conda/bin/conda", "run", "-n", "gadgetron", "--no-capture-output", "/opt/start_supervisor"]