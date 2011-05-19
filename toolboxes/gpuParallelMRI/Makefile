GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

DLLNAME=gpuparallelmri

LIBSOURCES=

HEADERFILES=\
	b1_map.h \
	htgrappa.h \
	cgOperatorSense.h \
	cgOperatorCartesianSense.h \
	cgOperatorNonCartesianSense.h \
	vector_td_identity_operator.h

CUDASOURCES=\
	b1_map.cu \
	cgOperatorSense.cu \
	cgOperatorCartesianSense.cu \
	cgOperatorNonCartesianSense.cu \
	htgrappa.cu \
	cposv_wrapper.cu \


LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

#LIBLDFLAGS += -lgpucore -lgpunfft -L../../dependencies/magma_1.0.0-rc4/lib -lmagma -lmagmablas -lf77blas -latlas -lcblas -lf2c -lcublas -lm -llapack 
#APPLDFLAGS += -lgpucore -lgpunfft -L../../dependencies/magma_1.0.0-rc4/lib -lmagma -lmagmablas -lf77blas -latlas -lcblas -lf2c -lcublas -lm -llapack 
LIBLDFLAGS += -lgpucore -lgpunfft -L../../dependencies/magma_1.0.0-rc4/lib -lmagma -lmagmablas
APPLDFLAGS += -lgpucore -lgpunfft -L../../dependencies/magma_1.0.0-rc4/lib -lmagma -lmagmablas
LIBLDFLAGS += -lgpucore -lgpunfft -L../../dependencies/magma_1.0.0-rc4/lib -latlas -lcblas -lcublas -lm -llapack 
APPLDFLAGS += -lgpucore -lgpunfft -L../../dependencies/magma_1.0.0-rc4/lib -latlas -lcblas -lcublas -lm -llapack 

CXXFLAGS += 
NVCCFLAGS += -I../../dependencies/magma_1.0.0-rc4/include
LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)

all: $(LIBFILE)

test_b1: $(LIBCUDAOBJECTS) b1map_test.cuo
	$(CXX) $(LIBOBJECTS) $(LIBCUDAOBJECTS) b1map_test.cuo $(CUDALIBRARIES) $(APPLDFLAGS) -o b1map_test 

test_htgrappa: $(LIBCUDAOBJECTS) htgrappa_test.o
	$(CXX) $(LIBOBJECTS) $(LIBCUDAOBJECTS) htgrappa_test.o $(CUDALIBRARIES) $(APPLDFLAGS) -o htgrappa_test 

$(LIBFILE): $(LIBOBJECTS) $(LIBCUDAOBJECTS) Makefile
	$(CXX) $(LIBLDFLAGS) $(LIBOBJECTS) $(LIBCUDAOBJECTS)  $(CUDALIBRARIES) -o $@
	mkdir -p $(GADGETRONHOME)/lib
	mkdir -p $(GADGETRONHOME)/inc
	cp $(LIBFILE) $(GADGETRONHOME)/lib/
	cp $(HEADERFILES) $(GADGETRONHOME)/inc/

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
