GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

DLLNAME=solvers

LIBHEADERS=\
	solver.h \
	cgSolver.h \
	cuCGSolver.h \
	sbSolver.h \
	cuSBSolver.h \
	cgPreconditioner.h \
	cuCGPreconditioner.h \
	cuCGPrecondWeights.h \
	matrixOperator.h \
	cuMatrixOperator_macros.h \
	identityOperator.h \
	cuIdentityOperator.h \
	imageOperator.h \
	cuImageOperator.h \
	partialDerivativeOperator.h \
	cuPartialDerivativeOperator.h \
	convolutionOperator.h \
	cuConvolutionOperator.h

LIBSOURCES=

CUDASOURCES=\
	cuCGPrecondWeights.cu \
	cuPartialDerivativeOperator.cu

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


