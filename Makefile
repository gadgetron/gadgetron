#Gadgetron Makefile

GADGETRONHOME=.

MACHINE   := $(shell uname -m)
KERNEL    := $(shell uname -s)

DLLEXTENSION=so
ifeq ($(KERNEL), Darwin)
DLLEXTENSION=dylib
endif


#CUDASDK=$(HOME)/NVIDIA_GPU_Computing_SDK
#CUDALIBRARIES=-L/usr/local/cuda/lib64 -L$(CUDASDK)/shared/lib -L$(CUDASDK)/C/lib -L$(CUDASDK)/C/common/lib/linux
#ifeq ($(KERNEL), Darwin)
#CUDALIBRARIES=-L/usr/local/cuda/lib -L$(CUDASDK)/shared/lib -L$(CUDASDK)/C/lib -L$(CUDASDK)/C/common/lib/darwin
#endif

HEADERS=\
	GadgetContainerMessage.h \
	Gadget.h \
	gadgetheaders.h \
	GadgetServerAcceptor.h \
	GadgetStreamConfiguratorFactory.h \
	GadgetStreamConfigurator.h \
	GadgetStreamController.h \
	Gadgetron.h

EXESOURCES=\
	main.cpp \
	GadgetStreamConfiguratorFactory.cpp \
	GadgetStreamController.cpp \
	GadgetServerAcceptor.cpp

LIBSOURCES=\
	GadgetStreamConfigurator.cpp 

EXEOBJECTS=$(EXESOURCES:.cpp=.o)
LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

CXX=g++
ifeq ($(KERNEL), Darwin)
CXX=g++ -m32 -arch i386
endif

EXELDFLAGS= -L. -L$(GADGETRONHOME)/lib -lACE -lgadgetron -lgadgetroncore -lgadgetrontools -lgadgetrongpucg 

LIBLDFLAGS= -shared -lACE 

CXXFLAGS=-c -fPIC -Wall -I. -I$(GADGETRONHOME)/inc/ -I./gadgettools/ -g #-DACE_NTRACE=0

EXECUTABLE=gadgetron
LIBFILE=libgadgetron.$(DLLEXTENSION)

all: $(LIBFILE) $(EXECUTABLE)

$(EXECUTABLE): $(EXEOBJECTS) Makefile
	$(CXX) $(EXELDFLAGS) $(EXEOBJECTS) -o $@

$(LIBFILE): $(LIBOBJECTS) Makefile
	$(CXX) $(LIBLDFLAGS) $(LIBOBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -rf *~
	rm -rf $(EXECUTABLE)
	rm -rf $(LIBFILE)
	rm -rf *.o
	rm -rf *.cplx
	rm -rf *.real

install:
	mkdir -p inc
	mkdir -p lib
	cp $(HEADERS) inc/
	cp $(LIBFILE) lib/
	cp $(EXECUTABLE) lib/
