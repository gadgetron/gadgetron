GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

DLLNAME=cuNFFT

HEADERFILES=\
	NFFT.h

CUDASOURCES=\
	NFFT.cu

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)
LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)
APPLDFLAGS += -lgpucore -lcuNFFT
LIBLDFLAGS += -lgpucore

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
	rm -rf main.cuo
	rm -rf test

cleandata:
	rm -rf *.cplx
	rm -rf *.real
