GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

DLLNAME=cucg

LIBHEADERS=\
	cuCG.h \
	cuCGMatrixOperator.h \
	cuCGPreconditioner.h

LIBSOURCES=

CUDASOURCES=\
	cuCG.cu \
	cuCGPreconditioner.cu

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)
LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)
LIBLDFLAGS +=

LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)

all: $(LIBFILE)

$(LIBFILE): $(LIBOBJECTS) $(LIBCUDAOBJECTS) Makefile
	$(CXX) $(LIBLDFLAGS) $(LIBOBJECTS) $(LIBCUDAOBJECTS)  $(CUDALIBRARIES) -o $@
	mkdir -p $(GADGETRONHOME)/lib
	mkdir -p $(GADGETRONHOME)/inc
	cp $(LIBFILE) $(GADGETRONHOME)/lib/
	cp $(LIBHEADERS) $(GADGETRONHOME)/inc/

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


