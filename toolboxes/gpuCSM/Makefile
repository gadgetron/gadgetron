GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

DLLNAME=cuCSM

LIBSOURCES=

HEADERFILES=\
	b1_map.hcu

CUDASOURCES=\
	b1_map.cu

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)


LIBLDFLAGS += -lcundarray -lcuNFFT
APPLDFLAGS += -lcundarray -lcuNFFT

LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)

all: $(LIBFILE)

test: $(LIBCUDAOBJECTS) main.cuo
	$(CXX) $(LIBOBJECTS) $(LIBCUDAOBJECTS) main.cuo $(CUDALIBRARIES) $(APPLDFLAGS) -o test 


$(LIBFILE): $(LIBOBJECTS) $(LIBCUDAOBJECTS) Makefile
	$(CXX) $(LIBLDFLAGS) $(LIBOBJECTS) $(LIBCUDAOBJECTS)  $(CUDALIBRARIES) -o $@
	mkdir -p $(GADGETRONHOME)/lib
	mkdir -p $(GADGETRONHOME)/inc
	cp $(LIBFILE) $(GADGETRONHOME)/lib/
	cp $(HEADERFILES) $(GADGETRONHOME)/inc/

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

%.cuo: %.cu
	$(NVCC) $(NVCCFLAGS) -o $@ $<

clean:
	rm -rf *~
	rm -rf $(LIBFILE)
	rm -rf $(LIBOBJECTS)
	rm -rf $(LIBCUDAOBJECTS)
	rm -rf test

cleandata:
	rm -rf *.cplx
	rm -rf *.real
