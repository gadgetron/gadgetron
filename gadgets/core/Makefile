#Gadgetron Makefile

UNAME := $(shell uname)

GADGETRONHOME=../..

DLLNAME=gadgetroncore

DLLEXTENSION=so

ifeq ($(UNAME), Darwin)
DLLEXTENSION=dylib
endif

CONFIG_FILES=\
	default.xml

LIBSOURCES=\
	AcquisitionPassthroughGadget.cpp \
	AcquisitionFinishGadget.cpp \
	AccumulatorGadget.cpp \
	FFTGadget.cpp \
	ImageFinishGadget.cpp \
	CropAndCombineGadget.cpp \
	ImageWriterGadget.cpp \
	MRIAcquisitionReader.cpp \
	MRIImageWriter.cpp 

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

CXX=g++

EXELDFLAGS= -lACE -lfftw3 -lfftw3f

CXXFLAGS=-c -fPIC -Wall -I. -I$(GADGETRONHOME)/inc -I$(GADGETRONHOME) -g #-DACE_NTRACE=0

ifeq ($(UNAME), Darwin)
CXXFLAGS=-m32 -arch i386 -c -fPIC -Wall -I. -I$(GADGETRONHOME)/inc -I$(GADGETRONHOME) -g #-DACE_NTRACE=0
endif 

LIBLDFLAGS= -shared -L$(GADGETRONHOME)/lib -lACE -lfftw3 -lfftw3f -ltinyxml
ifeq ($(UNAME), Darwin)
LIBLDFLAGS= -arch i386 -shared -L$(GADGETRONHOME)/lib -lACE -lfftw3 -lfftw3f -ltinyxml
endif

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
	mkdir -p $(GADGETRONHOME)/config
	cp $(LIBFILE) $(GADGETRONHOME)/lib/
	cp $(CONFIG_FILES) $(GADGETRONHOME)/config/
