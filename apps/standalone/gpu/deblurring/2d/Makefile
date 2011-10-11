GADGETRONHOME=../../../../..
include $(GADGETRONHOME)/Makefile.inc

APPLDFLAGS += -lgpucore -lsolvers -lhostutils

default: blur.o cg_deblur.o sb_deblur.o
	$(CXX) blur.o $(CUDALIBRARIES) $(APPLDFLAGS) -o blur
	$(CXX) cg_deblur.o $(CUDALIBRARIES) $(APPLDFLAGS) -o cg_deblur
	$(CXX) sb_deblur.o $(CUDALIBRARIES) $(APPLDFLAGS) -o sb_deblur

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

%.cuo: %.cu
	$(NVCC) $(NVCCFLAGS) -o $@ $<

clean:
	rm -rf *~
	rm -rf *.o
	rm blur
	rm cg_deblur
	rm sb_deblur
