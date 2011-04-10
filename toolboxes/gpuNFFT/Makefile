#Gadgetron Makefile

MACHINE   := $(shell uname -m)
KERNEL    := $(shell uname -s)


GADGETRONHOME=../..
CUDASDK=$(HOME)/NVIDIA_GPU_Computing_SDK
CUDAINCLUDES=\
	-I. \
	-I../../toolboxes/ndarray \
	-I../../toolboxes/gpundarray \
	-I/usr/local/cuda/include \
	-I$(CUDASDK)/C/common/inc \
	-I$(CUDASDK)/shared/inc

CUDAHOME=/usr/local/cuda

DLLNAME=cuNFFT

DLLEXTENSION=so

ifeq ($(KERNEL), Darwin)
DLLEXTENSION=dylib
endif

LIBSOURCES=

HEADERFILES=\
	NFFT.hcu

CUDASOURCES=\
	NFFT.cu \
	ndarray_device_utilities.cu

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

CXX=g++
ifeq ($(KERNEL), Darwin)
CXX=g++ -fPIC -arch i386 -m32 -Xlinker -rpath /usr/local/cuda/lib -O2
endif 

CXXFLAGS=-c -fPIC -Wall -I. \
	-I$(GADGETRONHOME)/inc \
	-I$(GADGETRONHOME)/gadgets/core \
	$(CUDAINCLUDES) \
	-DGADGETS_BUILD_DLL -O2

CXXFLAGS_MAIN=-c --Wall -I. \
	-I$(GADGETRONHOME)/inc \
	-I$(GADGETRONHOME)/gadgets/core \
	$(CUDAINCLUDES) \
	-O2

LIBLDFLAGS= -shared -lACE -L$(GADGETRONHOME)/lib -ltinyxml -lcundarray
MAINLDFLAGS= -lACE -L$(GADGETRONHOME)/lib -ltinyxml -lcundarray

NVCC=$(CUDAHOME)/bin/nvcc
GENCODE_SM10 := -gencode=arch=compute_10,code=\"sm_10,compute_10\"
GENCODE_SM20 := -gencode=arch=compute_20,code=\"sm_20,compute_20\"

GENCODE_SM := $(GENCODE_SM20)
ifeq ($(KERNEL), Darwin)
GENCODE_SM := $(GENCODE_SM10)
endif

NVCCFLAGS= $(GENCODE_SM) --compiler-options -fno-strict-aliasing $(CUDAINCLUDES) -O2 -c
ifeq ($(KERNEL), Darwin)
NVCCFLAGS += -m32
endif

CUDALIBRARIES=-L$(CUDAHOME)/lib64 -L$(CUDASDK)/shared/lib -L$(CUDASDK)/C/lib -L$(CUDASDK)/C/common/lib/linux -lshrutil_$(MACHINE) -lcutil_$(MACHINE) -lcudart -lcublas -lcufft

ifeq ($(KERNEL), Darwin)
CUDALIBRARIES=-L$(CUDAHOME)/lib -L$(CUDASDK)/shared/lib -L$(CUDASDK)/C/lib -L$(CUDASDK)/C/common/lib/darwin -lshrutil_$(MACHINE) -lcutil_$(MACHINE) -lcudart -lcublas -lcufft
endif 

LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)

all: $(LIBFILE)

test: $(LIBCUDAOBJECTS)
	$(NVCC) $(NVCCFLAGS) -o main.cuo main.cu
	$(CXX) $(LIBOBJECTS) $(LIBCUDAOBJECTS) main.cuo $(CUDALIBRARIES) $(MAINLDFLAGS) -o test 


$(LIBFILE): $(LIBOBJECTS) $(LIBCUDAOBJECTS) Makefile
	$(CXX) $(LIBLDFLAGS) $(LIBOBJECTS) $(LIBCUDAOBJECTS)  $(CUDALIBRARIES) -o $@

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

%.cuo: %.cu
	$(NVCC) $(NVCCFLAGS) -o $@ $<

clean:
	rm -rf *~
	rm -rf $(LIBFILE)
	rm -rf $(LIBOBJECTS)
	rm -rf $(LIBCUDAOBJECTS)
	rm -rf test

cleandata:
	rm -rf *.cplx
	rm -rf *.real

install:
	mkdir -p $(GADGETRONHOME)/lib
	mkdir -p $(GADGETRONHOME)/inc
	cp $(LIBFILE) $(GADGETRONHOME)/lib/
	cp $(HEADERFILES) $(GADGETRONHOME)/inc/
