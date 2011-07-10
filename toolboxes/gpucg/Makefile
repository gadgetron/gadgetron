GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

DLLNAME=gpucg

LIBHEADERS=\
	cuCG.h \
	cuCGPreconditioner.h \
	cuCGPrecondWeights.h \
	cuCGMatrixOperator.h \
	cuCGIdentityOperator.h \
	cuCGImageOperator.h

LIBSOURCES=

CUDASOURCES=\
	cuCG.cu \
	cuCGPrecondWeights.cu \
	cuCGImageOperator.cu

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)
LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)
LIBLDFLAGS += -lgpucore

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


