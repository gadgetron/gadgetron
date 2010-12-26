#Gadgetron Makefile

UNAME := $(shell uname)

DLLEXTENSION=so

ifeq ($(UNAME), Darwin)
DLLEXTENSION=dylib
endif

EXESOURCES=main.cpp GadgetStreamController.cpp GadgetServerAcceptor.cpp GadgetStreamConfigurator.cpp GadgetStreamConfiguratorFactory.cpp DefaultConfigurator.cpp AcquisitionPassthroughGadget.cpp gadgettools/ConfigParser.cpp AcquisitionFinishGadget.cpp
LIBSOURCES=

EXEOBJECTS=$(EXESOURCES:.cpp=.o)
LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

CXX=g++
#CXXFLAGS=-c -Wall -I.
EXELDFLAGS= -lACE -lfftw3 -lfftw3f

CXXFLAGS=-c -fPIC -Wall -I. -I./gadgettools/ #-DACE_NTRACE=0
LIBLDFLAGS= -shared -lACE 

LIBFILE=libgadgets.$(DLLEXTENSION)
EXECUTABLE=gadgetron

all: $(EXECUTABLE) $(LIBFILE)

$(EXECUTABLE): $(EXEOBJECTS) Makefile
	$(CXX) $(EXELDFLAGS) $(EXEOBJECTS) -o $@

$(LIBFILE): $(LIBOBJECTS) Makefile
	$(CXX) $(LIBLDFLAGS) $(LIBOBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -DGADGETS_BUILD_DLL $< -o $@

clean:
	rm -rf *~
	rm -rf $(EXECUTABLE)
	rm -rf $(LIBFILE)
	rm -rf *.o
