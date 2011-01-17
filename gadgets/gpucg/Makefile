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

DLLNAME=gadgetrongpucg

DLLEXTENSION=so

ifeq ($(KERNEL), Darwin)
DLLEXTENSION=dylib
endif

LIBHEADERS=\
	GPUCGConfigurator.h \
	GPUCGGadget.h 

LIBSOURCES=\
	GPUCGConfigurator.cpp \
	GPUCGGadget.cpp \
	GPUCGGoldenRadial.cpp \
	GPUCGGoldenSpiral.cpp \
	NSense.cpp \
	preprocess.cpp \
	ndarray.cpp \
	vec.cpp \
	vds.cpp

CUDASOURCES=\
	NFFT.cu \
	FFT.cu \
	image_utilities.cu \
	preprocess_.cu \
	preprocess_radial.cu \
	preprocess_sense.cu \
	radixSort/radixsort.cu \
	batchfft.cu \
	NSense_.cu

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

CXX=g++
ifeq ($(KERNEL), Darwin)
CXX=g++ -m32 -arch i386
endif 

CXXFLAGS=-c -fPIC -Wall -I. \
	-I$(GADGETRONHOME)/gadgettools \
	-I$(GADGETRONHOME)/inc \
	$(CUDAINCLUDES) \
	-g #-DACE_NTRACE=0

LIBLDFLAGS= -shared -lACE -L$(GADGETRONHOME)/lib -lgadgetron -lgadgetroncore -lgadgetrontools 

NVCC=$(CUDAHOME)/bin/nvcc
NVCCFLAGS= --compiler-options -fno-strict-aliasing -c $(CUDAINCLUDES) -Xcompiler -fPIC


CUDALIBRARIES=-L$(CUDAHOME)/lib64 -L$(CUDASDK)/shared/lib -L$(CUDASDK)/C/lib -L$(CUDASDK)/C/common/lib/linux -lshrutil_$(MACHINE) -lcutil_$(MACHINE)  -lcudpp_$(MACHINE) -lcudart -lcublas -lcufft

ifeq ($(KERNEL), Darwin)
CUDALIBRARIES=-L$(CUDAHOME)/lib -L$(CUDASDK)/shared/lib -L$(CUDASDK)/C/lib -L$(CUDASDK)/C/common/lib/darwin -lshrutil_$(MACHINE) -lcutil_$(MACHINE)  -lcudpp_$(MACHINE) -lcudart -lcublas -lcufft
endif 

LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)

all: $(LIBFILE)

test:
	echo $(LIBCUDAOBJECTS)


$(LIBFILE): $(LIBOBJECTS) $(LIBCUDAOBJECTS) Makefile
	$(CXX) $(LIBLDFLAGS) $(LIBOBJECTS) $(LIBCUDAOBJECTS)  $(CUDALIBRARIES) -o $@

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) -DGADGETS_BUILD_DLL $< -o $@

%.cuo: %.cu
	$(NVCC) $(NVCCFLAGS) -DGADGETS_BUILD_DLL -o $@ $<

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
	cp $(LIBHEADERS) $(GADGETRONHOME)/inc/
	cp $(LIBFILE) $(GADGETRONHOME)/lib/
