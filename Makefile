#Gadgetron Makefile

GADGETRONHOME=.



MACHINE   := $(shell uname -m)
KERNEL    := $(shell uname -s)


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
	GadgetStreamController.cpp \
	GadgetServerAcceptor.cpp \
	GadgetStreamConfigurator.cpp \
	GadgetStreamConfiguratorFactory.cpp 

EXEOBJECTS=$(EXESOURCES:.cpp=.o)

CXX=g++
EXELDFLAGS= -L$(GADGETRONHOME)/lib -lACE -lgadgetroncore -lgadgetrontools -lgadgetrongpucg 

CXXFLAGS=-c -fPIC -Wall -I. -I$(GADGETRONHOME)/inc/ -I./gadgettools/ -g #-DACE_NTRACE=0

EXECUTABLE=gadgetron

all: $(EXECUTABLE)

$(EXECUTABLE): $(EXEOBJECTS) Makefile
	$(CXX) $(EXELDFLAGS) $(EXEOBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -rf *~
	rm -rf $(EXECUTABLE)
	rm -rf *.o
	rm -rf *.cplx
	rm -rf *.real

install:
	mkdir -p inc
	mkdir -p lib
	cp $(HEADERS) inc/
	cp $(EXECUTABLE) lib/