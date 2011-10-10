GADGETRONHOME=../../../../../../../..
include $(GADGETRONHOME)/Makefile.inc

EXECUTABLE=radial_ktsense

HEADERFILES=

APPSOURCES=\
	main.cpp

CUDASOURCES=

APPOBJECTS=$(APPSOURCES:.cpp=.o)
CUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

APPLDFLAGS += -lhostutils -lgpucore -lgpunfft -lsolvers -lgpuparallelmri

test: main.o $(LIBCUDAOBJECTS)
	$(CXX) $(APPOBJECTS) $(LIBCUDAOBJECTS)  $(CUDALIBRARIES) $(APPLDFLAGS) -o $(EXECUTABLE) 
	mkdir -p $(GADGETRONHOME)/bin
	cp $(EXECUTABLE) $(GADGETRONHOME)/bin/

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

%.cuo: %.cu
	$(NVCC) $(NVCCFLAGS)  -o $@ $<

clean:
	rm -rf *~
	rm -rf $(EXECUTABLE)
	rm -rf $(APPOBJECTS)
	rm -rf $(CUDAOBJECTS)
