#Gadgetron Makefile

MACHINE   := $(shell uname -m)
KERNEL    := $(shell uname -s)


GADGETRONHOME=../..
CUDASDK=$(HOME)/NVIDIA_GPU_Computing_SDK
CUDAINCLUDES=\
	-I. \
	-I/usr/local/cuda/include \
	-I$(CUDASDK)/C/common/inc \
	-I$(CUDASDK)/shared/inc

CUDAHOME=/usr/local/cuda

DLLNAME=cundarray

DLLEXTENSION=so

ifeq ($(KERNEL), Darwin)
DLLEXTENSION=dylib
endif

LIBHEADERS=\
	cuNDArray.h \
	cuNDFFT.h

LIBSOURCES=\
	cuNDFFT.cpp \
	cuNDArray.cpp

CUDASOURCES=\
	cuNDArray_kernels.cu \
	FFT.cu


LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

CXX=g++
ifeq ($(KERNEL), Darwin)
CXX=g++ -m32 -arch i386
endif 

CXXFLAGS=-c -fPIC -Wall -I. \
	-I$(GADGETRONHOME)/inc \
	-I$(GADGETRONHOME)/gadgets/core \
	$(CUDAINCLUDES) \
	-DGADGETS_BUILD_DLL -g #-DACE_NTRACE=0

LIBLDFLAGS= -shared -lACE -L$(GADGETRONHOME)/lib -ltinyxml

NVCC=$(CUDAHOME)/bin/nvcc
ifeq ($(KERNEL), Darwin)
NVCC=$(CUDAHOME)/bin/nvcc -m32 -gencode=arch=compute_10,code=\"sm_10,compute_10\"
endif 

NVCCFLAGS= --compiler-options -fno-strict-aliasing -c $(CUDAINCLUDES) -Xcompiler -fPIC -DGADGETS_BUILD_DLL -I$(GADGETRONHOME)/inc

CUDALIBRARIES=-L$(CUDAHOME)/lib64 -L$(CUDASDK)/shared/lib -L$(CUDASDK)/C/lib -L$(CUDASDK)/C/common/lib/linux -lshrutil_$(MACHINE) -lcutil_$(MACHINE)  -lcudart -lcublas -lcufft

ifeq ($(KERNEL), Darwin)
CUDALIBRARIES=-L$(CUDAHOME)/lib -L$(CUDASDK)/shared/lib -L$(CUDASDK)/C/lib -L$(CUDASDK)/C/common/lib/darwin -lshrutil_$(MACHINE) -lcutil_$(MACHINE)  -lcudart -lcublas -lcufft
endif 

LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)

all: $(LIBFILE)

$(LIBFILE): $(LIBOBJECTS) $(LIBCUDAOBJECTS) Makefile
	$(CXX) $(LIBLDFLAGS) $(LIBOBJECTS) $(LIBCUDAOBJECTS)  $(CUDALIBRARIES) -o $@

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

%.cuo: %.cu
	$(NVCC) $(NVCCFLAGS)  -o $@ $<

clean:
	rm -rf *~
	rm -rf $(LIBFILE)
	rm -rf $(LIBOBJECTS)
	rm -rf $(LIBCUDAOBJECTS)
	rm -rf *.cplx
	rm -rf *.real

install:
	mkdir -p $(GADGETRONHOME)/lib
	mkdir -p $(GADGETRONHOME)/inc
	cp $(LIBFILE) $(GADGETRONHOME)/lib/
	cp $(LIBHEADERS) $(GADGETRONHOME)/inc/

