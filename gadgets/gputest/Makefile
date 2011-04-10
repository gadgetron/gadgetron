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

DLLNAME=gadgetrongputest

DLLEXTENSION=so

ifeq ($(KERNEL), Darwin)
DLLEXTENSION=dylib
endif

CONFIG_FILES=\
	gputest.xml

LIBSOURCES=
#LIBSOURCES=\
	GPUCGGadget.cpp \
	GPUCGGoldenRadial.cpp \
	GPUCGGoldenSpiral.cpp \
	NSense.cpp \
	preprocess.cpp \
	ndarray.cpp \
	vec.cpp \
	vds.cpp

CUDASOURCES=\
	cgOperatorSense.cu \
	cgOperatorCartesianSense.cu \
	cgOperatorNonCartesianSense.cu


LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

CXX=g++ -g
ifeq ($(KERNEL), Darwin)
CXX=g++ -m32 -arch i386
endif 

CXXFLAGS=-c -fPIC -Wall -I. \
	-I$(GADGETRONHOME)/inc \
	-I$(GADGETRONHOME)/gadgets/core \
	$(CUDAINCLUDES) \
	-DGADGETS_BUILD_DLL -g #-DACE_NTRACE=0

CXXFLAGS_MAIN=-c -Wall -I. \
	-I$(GADGETRONHOME)/inc \
	-I$(GADGETRONHOME)/gadgets/core \
	$(CUDAINCLUDES)

LIBLDFLAGS= -shared -lACE -L$(GADGETRONHOME)/lib -ltinyxml
MAINLDFLAGS= -lACE -L$(GADGETRONHOME)/lib -ltinyxml -lcundarray -lcucg -lcuNFFT

NVCC=$(CUDAHOME)/bin/nvcc
NVCCFLAGS= -g --compiler-options -fno-strict-aliasing -c $(CUDAINCLUDES) -Xcompiler -fPIC -DGADGETS_BUILD_DLL -I$(GADGETRONHOME)/inc


CUDALIBRARIES=-L$(CUDAHOME)/lib64 -L$(CUDASDK)/shared/lib -L$(CUDASDK)/C/lib -L$(CUDASDK)/C/common/lib/linux -lshrutil_$(MACHINE) -lcutil_$(MACHINE)  -lcudart -lcublas -lcufft

ifeq ($(KERNEL), Darwin)
CUDALIBRARIES=-L$(CUDAHOME)/lib -L$(CUDASDK)/shared/lib -L$(CUDASDK)/C/lib -L$(CUDASDK)/C/common/lib/darwin -lshrutil_$(MACHINE) -lcutil_$(MACHINE)  -lcudart -lcublas -lcufft
endif 

LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)

all: $(LIBFILE)

test: main.o $(LIBCUDAOBJECTS)
	$(CXX) $(CXXFLAGS_MAIN) main.cpp
	$(CXX) main.o $(LIBOBJECTS) $(LIBCUDAOBJECTS)  $(CUDALIBRARIES) $(MAINLDFLAGS) -o test 


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
#	rm -rf *.cplx
#	rm -rf *.real
	rm -rf main.o 

install:
	mkdir -p $(GADGETRONHOME)/lib
	mkdir -p $(GADGETRONHOME)/config
	cp $(LIBFILE) $(GADGETRONHOME)/lib/
	cp $(CONFIG_FILES) $(GADGETRONHOME)/config/

