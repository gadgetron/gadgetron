
default: build

configure:
    mkdir -p build/ ;\
    cd build/ ;\
    cmake -GNinja -D CMAKE_BUILD_TYPE=Release \
        -D USE_CUDA=ON -D USE_MKL=ON \
        -D BUILD_DOCUMENTATION=OFF \
        -D CMAKE_INSTALL_PREFIX=${CONDA_PREFIX} \
        ../

build: configure
    cd build && ninja

install: build
    cd build && ninja install

test: unit-test e2e-test

unit-test: build
    cd build && ctest

e2e-test: install
    cd test/e2e && pytest
