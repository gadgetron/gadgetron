GADGETRONHOME=../../../../../../../..
include $(GADGETRONHOME)/Makefile.inc

EXECUTABLE_CG=radial_sense_cg
EXECUTABLE_SB=radial_sense_sb

APPLDFLAGS += -lhostutils -lgpucore -lgpunfft -lsolvers -lgpuparallelmri

default: main_cg.o main_sb.o
	$(CXX) main_cg.o $(CUDALIBRARIES) $(APPLDFLAGS) -o $(EXECUTABLE_CG) 
	$(CXX) main_sb.o $(CUDALIBRARIES) $(APPLDFLAGS) -o $(EXECUTABLE_SB) 
	mkdir -p $(GADGETRONHOME)/bin
	cp $(EXECUTABLE_CG) $(GADGETRONHOME)/bin/
	cp $(EXECUTABLE_SB) $(GADGETRONHOME)/bin/

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

%.cuo: %.cu
	$(NVCC) $(NVCCFLAGS)  -o $@ $<

clean:
	rm -rf *~
	rm -rf *.o
	rm -rf $(EXECUTABLE_SB)
	rm -rf $(EXECUTABLE_CG)
	rm -rf $(APPOBJECTS)
	rm -rf $(CUDAOBJECTS)
