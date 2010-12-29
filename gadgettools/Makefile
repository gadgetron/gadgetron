#Gadgetron Makefile

UNAME := $(shell uname)

GADGETRONHOME=..

DLLNAME=gadgetrontools

DLLEXTENSION=so

ifeq ($(UNAME), Darwin)
DLLEXTENSION=dylib
endif

LIBHEADERS=\
	ConfigParser.h \
	FFT.h \
	GadgetSocketReceiver.h \
	GadgetSocketSender.h \
	NDArray.h

LIBSOURCES=\
	ConfigParser.cpp

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

CXX=g++
EXELDFLAGS= -lACE

CXXFLAGS=-c -fPIC -Wall -I. -I$(GADGETRONHOME) -g #-DACE_NTRACE=0
LIBLDFLAGS= -shared -lACE 

LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)

all: $(LIBFILE)

$(LIBFILE): $(LIBOBJECTS) Makefile
	$(CXX) $(LIBLDFLAGS) $(LIBOBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -DGADGETS_BUILD_DLL $< -o $@

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
