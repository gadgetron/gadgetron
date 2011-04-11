GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

LIBHEADERS=\
	cuNDArray.h \
	cuNDFFT.h

LIBSOURCES=\
	cuNDFFT.cpp \
	cuNDArray.cpp

CUDASOURCES=\
	cuNDArray_kernels.cu 

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)
LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)
DLLNAME=cundarray
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


