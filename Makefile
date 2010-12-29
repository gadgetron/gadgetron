#Gadgetron Makefile

GADGETRONHOME=.

UNAME := $(shell uname)

HEADERS=\
	GadgetContainerMessage.h \
	Gadget.h \
	gadgetheaders.h \
	GadgetServerAcceptor.h \
	GadgetStreamConfiguratorFactory.h \
	GadgetStreamConfigurator.h \
	GadgetStreamController.h

EXESOURCES=\
	main.cpp \
	GadgetStreamController.cpp \
	GadgetServerAcceptor.cpp \
	GadgetStreamConfigurator.cpp \
	GadgetStreamConfiguratorFactory.cpp 

EXEOBJECTS=$(EXESOURCES:.cpp=.o)

CXX=g++
EXELDFLAGS= -L$(GADGETRONHOME)/lib -lACE -lfftw3 -lfftw3f -lgadgetroncore -lgadgetrontools

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
	cp $(HEADERS) inc/
	cp $(EXECUTABLE) lib/