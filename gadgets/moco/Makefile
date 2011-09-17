GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

CONFIG_FILES=\
	moco.xml

LIBSOURCES=\

CUDASOURCES=

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)
LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

DLLNAME=gadgetronmoco

LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)
LIBLDFLAGS += -lgpuparallelmri -lfftw3 -lfftw3f -lgpucore -lACE -ltinyxml

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
	rm -rf *.cplx
	rm -rf *.real
