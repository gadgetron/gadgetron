GADGETRONHOME=../../..
include $(GADGETRONHOME)/Makefile.inc

EXECUTABLE=gputestapp

APPSOURCES=\
	main.cpp
CUDASOURCES=

APPOBJECTS=$(APPSOURCES:.cpp=.o)
CUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

APPLDFLAGS += -lgpucore -lgpunfft -lgpucg -lgpuparallelmri
APPLDFLAGS += -L$(GADGETRONHOME)/dependencies/magma_1.0.0-rc4/lib -lmagma -lmagmablas

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
#	rm -rf *.cplx
#	rm -rf *.real


