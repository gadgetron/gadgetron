#Gadgetron Makefile

GADGETRONHOME=.

MACHINE   := $(shell uname -m)
KERNEL    := $(shell uname -s)

DLLEXTENSION=so
ifeq ($(KERNEL), Darwin)
DLLEXTENSION=dylib
endif

HEADERS=\
	GadgetContainerMessage.h \
	Gadget.h \
	GadgetServerAcceptor.h \
	GadgetStreamController.h \
	Gadgetron.h \
	GadgetronExport.h \
	GadgetMessageInterface.h \

EXESOURCES=\
	main.cpp \
	GadgetStreamController.cpp \
	GadgetServerAcceptor.cpp

EXEOBJECTS=$(EXESOURCES:.cpp=.o)
LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

CXX=g++
ifeq ($(KERNEL), Darwin)
CXX=g++ -m32 -arch i386
endif

EXELDFLAGS=-lACE -ltinyxml

LIBLDFLAGS= -shared -lACE 

CXXFLAGS=-c -fPIC -Wall -I.  -I$(GADGETRONHOME)/toolboxes/gadgettools -g #-DACE_NTRACE=0

EXECUTABLE=gadgetron
LIBFILE=libgadgetron.$(DLLEXTENSION)

all: $(EXECUTABLE)

$(EXECUTABLE): $(EXEOBJECTS) Makefile
	$(CXX) $(EXELDFLAGS) $(EXEOBJECTS) -o $@

#$(LIBFILE): $(LIBOBJECTS) Makefile
#	$(CXX) $(LIBLDFLAGS) $(LIBOBJECTS) -o $@

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
	cp $(EXECUTABLE) lib/
