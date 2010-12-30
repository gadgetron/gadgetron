#Gadgetron Makefile

UNAME := $(shell uname)

GADGETRONHOME=../..

DLLNAME=gadgetroncore

DLLEXTENSION=so

ifeq ($(UNAME), Darwin)
DLLEXTENSION=dylib
endif

LIBHEADERS=\
	DefaultConfigurator.h \
	AcquisitionPassthroughGadget.h \
	AcquisitionFinishGadget.h \
	AccumulatorGadget.h \
	FFTGadget.h \
	ImageFinishGadget.h \
	CropAndCombineGadget.h \
	ImageWriterGadget.h

LIBSOURCES=\
	DefaultConfigurator.cpp \
	AcquisitionPassthroughGadget.cpp \
	AcquisitionFinishGadget.cpp \
	AccumulatorGadget.cpp \
	FFTGadget.cpp \
	ImageFinishGadget.cpp \
	CropAndCombineGadget.cpp \
	ImageWriterGadget.cpp

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

CXX=g++
EXELDFLAGS= -lACE -lfftw3 -lfftw3f

CXXFLAGS=-c -fPIC -Wall -I. -I$(GADGETRONHOME)/gadgettools -I$(GADGETRONHOME) -g #-DACE_NTRACE=0
LIBLDFLAGS= -shared -lACE -lfftw3 -lfftw3f

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
