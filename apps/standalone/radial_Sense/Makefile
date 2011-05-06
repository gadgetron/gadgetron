GADGETRONHOME=../../..
include $(GADGETRONHOME)/Makefile.inc

EXECUTABLE=radial_Sense

HEADERFILES=\
        cgOperatorSense.h \
	cgOperatorNonCartesianSense.h

APPSOURCES=\
	main.cpp

CUDASOURCES=\
	cgOperatorSense.cu \
	cgOperatorNonCartesianSense.cu

APPOBJECTS=$(APPSOURCES:.cpp=.o)
CUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

APPLDFLAGS += -lhostutils -lgpucore -lgpunfft -lgpucsm -lgpucg

test: main.o $(LIBCUDAOBJECTS)
	$(CXX) $(APPOBJECTS) $(LIBCUDAOBJECTS)  $(CUDALIBRARIES) $(APPLDFLAGS) -o $(EXECUTABLE) 
	mkdir -p $(GADGETRONHOME)/bin
	cp $(EXECUTABLE) $(GADGETRONHOME)/bin/
	cp $(HEADERFILES) $(GADGETRONHOME)/inc/

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

%.cuo: %.cu
	$(NVCC) $(NVCCFLAGS)  -o $@ $<

clean:
	rm -rf *~
	rm -rf $(EXECUTABLE)
	rm -rf $(APPOBJECTS)
	rm -rf $(CUDAOBJECTS)
