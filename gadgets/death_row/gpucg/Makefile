GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

CONFIG_FILES=\
	gpucg_radial_single.xml \
	gpucg_radial_3_slices.xml \
	gpucg_spiral_single.xml \
	gpucg_spiral_3_slices.xml

LIBSOURCES=\
	GPUCGGadget.cpp \
	GPUCGGoldenRadial.cpp \
	GPUCGGoldenSpiral.cpp \
	NSense.cpp \
	preprocess.cpp \
	ndarray.cpp \
	vec.cpp \
	vds.cpp

CUDASOURCES=\
	NFFT.cu \
	FFT.cu \
	image_utilities.cu \
	preprocess_.cu \
	preprocess_radial.cu \
	preprocess_sense.cu \
	radixSort/radixsort.cu \
	batchfft.cu \
	NSense_.cu

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)
LIBCUDAOBJECTS=$(CUDASOURCES:.cu=.cuo)

DLLNAME=gadgetrongpucg
LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)

CXXFLAGS += -I../core

all: $(LIBFILE)

$(LIBFILE): $(LIBOBJECTS) $(LIBCUDAOBJECTS) Makefile
	$(CXX) $(LIBLDFLAGS) $(LIBOBJECTS) $(LIBCUDAOBJECTS)  $(CUDALIBRARIES) -o $@
	mkdir -p $(GADGETRONHOME)/lib
	mkdir -p $(GADGETRONHOME)/config
	cp $(LIBFILE) $(GADGETRONHOME)/lib/
	cp $(CONFIG_FILES) $(GADGETRONHOME)/config/

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

