GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

CONFIG_FILES=\
	grappa.xml

LIBSOURCES=\
	GrappaGadget.cpp \
	GrappaCalibrationBuffer.cpp \
	GrappaWeights.cpp \
	GrappaWeightsCalculator.cpp

CUDASOURCES=

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)
LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

DLLNAME=gadgetrongrappa
LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)

CXXFLAGS += -I../core

all: $(LIBFILE)

$(LIBFILE): $(LIBOBJECTS) $(LIBCUDAOBJECTS) Makefile
	$(CXX) $(LIBLDFLAGS) $(LIBOBJECTS) $(LIBCUDAOBJECTS)  $(CUDALIBRARIES) -o $@
	mkdir -p $(GADGETRONHOME)/lib
	mkdir -p $(GADGETRONHOME)/config
	cp $(LIBFILE) $(GADGETRONHOME)/lib/
	cp $(CONFIG_FILES) $(GADGETRONHOME)/config/

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

%.cuo: %.cu
	$(NVCC) $(NVCCFLAGS)  -o $@ $<

clean:
	rm -rf *~
	rm -rf $(LIBFILE)
	rm -rf $(LIBOBJECTS)
	rm -rf $(LIBCUDAOBJECTS)
	rm -rf *.cplx
	rm -rf *.real

