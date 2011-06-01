MACHINE   := $(shell uname -m)
KERNEL    := $(shell uname -s)
CUDAHOME=/usr/local/cuda

CUDASDK=$(HOME)/NVIDIA_GPU_Computing_SDK
CUDAINCLUDES=\
	-I/usr/local/cuda/include \
	-I$(CUDASDK)/C/common/inc \
	-I$(CUDASDK)/shared/inc \
	-I/usr/local/cula/include 

DLLEXTENSION=so

ifeq ($(KERNEL), Darwin)
DLLEXTENSION=dylib
endif

LIBSOURCES=

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

CXX=g++ -g
ifeq ($(KERNEL), Darwin)
CXX=g++ -g -m32 -arch i386
endif 

NVCC=$(CUDAHOME)/bin/nvcc
GENCODE_SM10=-gencode=arch=compute_10,code=\"sm_10,compute_10\"
GENCODE_SM20=-gencode=arch=compute_20,code=\"sm_20,compute_20\"

GENCODE_SM=$(GENCODE_SM20)
ifeq ($(KERNEL), Darwin)
GENCODE_SM=$(GENCODE_SM10)
endif

NVCCFLAGS= $(GENCODE_SM) -g --compiler-options -fno-strict-aliasing -Xcompiler -fPIC -I$(GADGETRONHOME)/inc $(CUDAINCLUDES) -c
ifeq ($(KERNEL), Darwin)
NVCCFLAGS += -m32
endif

CXXFLAGS=-c -fPIC -Wall -I. \
	-I$(GADGETRONHOME)/inc \
	$(CUDAINCLUDES)


CUDALIBRARIES=-L$(CUDAHOME)/lib64 -L$(CUDASDK)/shared/lib -L$(CUDASDK)/C/lib -L$(CUDASDK)/C/common/lib/linux -L/usr/local/cula/lib64 -lcula -lshrutil_$(MACHINE) -lcutil_$(MACHINE)  -lcudart -lcublas -lcufft

ifeq ($(KERNEL), Darwin)
CUDALIBRARIES=-L$(CUDAHOME)/lib -L$(CUDASDK)/shared/lib -L$(CUDASDK)/C/lib -L$(CUDASDK)/C/common/lib/darwin -L/usr/local/cula/lib -lcula -lshrutil_$(MACHINE) -lcutil_$(MACHINE)  -lcudart -lcublas -lcufft
endif 

LIBLDFLAGS=-shared -L$(GADGETRONHOME)/lib
APPLDFLAGS=-L$(GADGETRONHOME)/lib
