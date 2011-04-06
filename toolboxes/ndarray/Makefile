#Gadgetron Makefile

UNAME := $(shell uname)

GADGETRONHOME=../..

DLLNAME=ndarray

DLLEXTENSION=so

ifeq ($(UNAME), Darwin)
DLLEXTENSION=dylib
endif

LIBHEADERS=\
	NDArray.h \
	hoNDArray.h \
	hoNDArray_fileio.h \
	FFT.h

#LIBSOURCES=\
	ConfigParser.cpp

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

EXELDFLAGS= -lACE

CXXFLAGS=-c -fPIC -Wall -I. -I$(GADGETRONHOME) -g #-DACE_NTRACE=0

ifeq ($(UNAME), Darwin)
CXXFLAGS=-arch i386 -c -fPIC -Wall -I. -I$(GADGETRONHOME) -g #-DACE_NTRACE=0
endif

LIBLDFLAGS= -shared -lACE 
ifeq ($(UNAME), Darwin)
LIBLDFLAGS= -arch i386 -shared -lACE 
endif

LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)

#all: $(LIBFILE)

install:
	mkdir -p $(GADGETRONHOME)/lib
	mkdir -p $(GADGETRONHOME)/inc
	cp $(LIBHEADERS) $(GADGETRONHOME)/inc/

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

