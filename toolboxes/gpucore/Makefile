GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

LIBHEADERS=\
	GPUTimer.h \
	check_CUDA.h \
	vector_td.h \
	vector_td_operators.h \
	vector_td_utilities.h \
	ndarray_vector_td_utilities.h \
	cuNDArray.h \
	cuNDFFT.h \
	radial_utilities.h \
	real_utilities.h \
	real_utilities_device.h

LIBSOURCES=\
	cuNDArray.cpp \
	cuNDFFT.cpp

CUDASOURCES= \
	ndarray_vector_td_utilities.cu \
	cuNDArray_kernels.cu \
	radial_utilities.cu

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)
LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)
DLLNAME=gpucore
LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)
LIBLDFLAGS +=

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


