GADGETRONHOME=../../../../../..
include $(GADGETRONHOME)/Makefile.inc

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)
LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

APPLDFLAGS += -lgpucore -lgpunfft -lhostutils

default: nfft_main.o nffth_main.o
	$(CXX) nfft_main.o $(CUDALIBRARIES) $(APPLDFLAGS) -o nfft 
	$(CXX) nffth_main.o $(CUDALIBRARIES) $(APPLDFLAGS) -o nffth

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

%.cuo: %.cu
	$(NVCC) $(NVCCFLAGS) -o $@ $<

clean:
	rm -rf *~
	rm -rf *.cuo
	rm -rf *.o
	rm -rf nfft
	rm -rf nffth
