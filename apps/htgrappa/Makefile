GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

EXECUTABLE=htgrappa

APPSOURCES=\
	main.cpp

CUDASOURCES=\
	htgrappa.cu \
	cposv_wrapper.cu

APPOBJECTS=$(APPSOURCES:.cpp=.o)
CUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

APPLDFLAGS += -lcuNFFT -lcundarray -lcucg -lgpucore -lcuCSM -L../../dependencies/magma_1.0.0-rc4/lib -lmagma -lmagmablas -lf77blas -latlas -lcblas -lf2c -lcublas -lm -llapack 

NVCCFLAGS += -I../../dependencies/magma_1.0.0-rc4/include

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


