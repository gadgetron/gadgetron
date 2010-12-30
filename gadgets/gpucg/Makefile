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
	NSense.cpp \
	preprocess.cpp \
	ndarray.cpp \
	vec.cpp

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

#LIBCUDAOBJECTS=$OBJS +=  $(patsubst %.cu,$(OBJDIR)/%.cu.o,$(notdir $(CUFILES)))(CUDASOURCES:.cu=.cu.o)

LIBCUDAOBJECTS=$(patsubst %.cu,%.cu.o,$(CUDASOURCES))

CXX=g++

CXXFLAGS=-c -fPIC -Wall -I. \
	-I$(GADGETRONHOME)/gadgettools \
	-I$(GADGETRONHOME)/inc \
	$(CUDAINCLUDES) \
	-g #-DACE_NTRACE=0

LIBLDFLAGS= -shared -lACE 

NVCC=$(CUDAHOME)/bin/nvcc
NVCCFLAGS=-c $(CUDAINCLUDES) -Xcompiler -fPIC 


LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)

all: $(LIBFILE)

$(LIBFILE): $(LIBOBJECTS) $(LIBCUDAOBJECTS) Makefile
	$(CXX) -m64 $(LIBLDFLAGS) $(LIBOBJECTS) $(LIBCUDAOBJECTS) -L$(CUDASDK)/shared/lib -L$(CUDASDK)/C/lib -L$(CUDASDK)/C/common/lib/linux -lshrutil_$(MACHINE) -lcutil_$(MACHINE) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -DGADGETS_BUILD_DLL $< -o $@

%.cu.o: %.cu
	$(NVCC) $(NVCCFLAGS) -DGADGETS_BUILD_DLL -o $@ $<

clean:
	rm -rf *~
	rm -rf $(LIBFILE)
	rm -rf *.o
	rm -rf *.cplx
	rm -rf *.real

install:
	mkdir -p $(GADGETRONHOME)/lib
	mkdir -p $(GADGETRONHOME)/inc
	cp $(LIBFILE) $(GADGETRONHOME)/lib/
	cp $(LIBHEADERS) $(GADGETRONHOME)/inc/
