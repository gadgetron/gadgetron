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

CONFIG_FILES=\
	gpucg_radial_single.xml \
	gpucg_radial_3_slices.xml \
	gpucg_spiral_single.xml \
	gpucg_spiral_3_slices.xml

LIBHEADERS=\
	GPUCGGadget.h 

LIBSOURCES=\
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
	-DGADGETS_BUILD_DLL -g #-DACE_NTRACE=0

LIBLDFLAGS= -shared -lACE -L$(GADGETRONHOME)/lib -lticpp

NVCC=$(CUDAHOME)/bin/nvcc
NVCCFLAGS= --compiler-options -fno-strict-aliasing -c $(CUDAINCLUDES) -Xcompiler -fPIC -DGADGETS_BUILD_DLL


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
	mkdir -p $(GADGETRONHOME)/config
	cp $(LIBHEADERS) $(GADGETRONHOME)/inc/
	cp $(LIBFILE) $(GADGETRONHOME)/lib/
	cp $(CONFIG_FILES) $(GADGETRONHOME)/config/

