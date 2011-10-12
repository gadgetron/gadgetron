GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

CONFIG_FILES=\
	radial_single.xml

LIBSOURCES=\
	GPUCGGadget.cpp \
	GPUCGGoldenRadial.cpp

CUDASOURCES=

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)
LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

DLLNAME=gadgetroncgsense

LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)
LIBLDFLAGS += -lgpuparallelmri -lgpucore -lgpunfft -lsolvers -lACE -ltinyxml

CXXFLAGS += -I../core

all: $(LIBFILE)

$(LIBFILE): $(LIBOBJECTS) $(LIBCUDAOBJECTS) $(CONFIG_FILES) Makefile
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
