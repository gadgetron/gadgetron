GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

DLLNAME=cuNFFT

HEADERFILES=\
	NFFT.hcu \
	ndarray_device_utilities.hcu \
	vector_utilities.hcu \
	floatd_operators.hcu \
	floatd_utilities.hcu \
	uintd_operators.hcu \
	uintd_utilities.hcu \
	intd_operators.hcu \
	intd_utilities.hcu \
	check_CUDA.h

CUDASOURCES=\
	NFFT.cu \
	ndarray_device_utilities.cu

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)
LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)
APPLDFLAGS += -lcuNFFT -lcundarray

all: $(LIBFILE)
	mkdir -p $(GADGETRONHOME)/lib
	mkdir -p $(GADGETRONHOME)/inc
	cp $(LIBFILE) $(GADGETRONHOME)/lib/
	cp $(HEADERFILES) $(GADGETRONHOME)/inc/

test: $(LIBCUDAOBJECTS) main.cuo
#	$(NVCC) $(NVCCFLAGS) -o main.cuo main.cu
	$(CXX) $(LIBOBJECTS) $(LIBCUDAOBJECTS) main.cuo $(CUDALIBRARIES) $(APPLDFLAGS) -o test 


$(LIBFILE): $(LIBOBJECTS) $(LIBCUDAOBJECTS) Makefile
	$(CXX) $(LIBLDFLAGS) $(LIBOBJECTS) $(LIBCUDAOBJECTS)  $(CUDALIBRARIES) -o $@

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

