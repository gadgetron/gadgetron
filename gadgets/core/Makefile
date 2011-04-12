GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

DLLNAME=gadgetroncore

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

LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)

LIBLDFLAGS += -lfftw3 -lfftw3f

all: $(LIBFILE)

$(LIBFILE): $(LIBOBJECTS) Makefile
	$(CXX) $(LIBLDFLAGS) $(LIBOBJECTS) -o $@
	mkdir -p $(GADGETRONHOME)/lib
	mkdir -p $(GADGETRONHOME)/config
	cp $(LIBFILE) $(GADGETRONHOME)/lib/
	cp $(CONFIG_FILES) $(GADGETRONHOME)/config/

.cpp.o:
	$(CXX) $(CXXFLAGS) -DGADGETS_BUILD_DLL $< -o $@

clean:
	rm -rf *~
	rm -rf $(LIBFILE)
	rm -rf *.o
	rm -rf *.cplx
	rm -rf *.real

