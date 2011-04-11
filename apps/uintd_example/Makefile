GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

CUDASOURCES=\
	uintd.cu \
	cuNDArray_new_permute.cu

LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)


APPLDFLAGS += -lcundarray -lcucg -lcuNFFT

test: main.o $(LIBCUDAOBJECTS)
	$(CXX) main.o $(LIBOBJECTS) $(LIBCUDAOBJECTS)  $(CUDALIBRARIES) $(APPLDFLAGS) -o test 


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
	rm -rf main.o 
	rm -rf test