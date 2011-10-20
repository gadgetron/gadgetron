GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

DLLNAME=gadgetronpython

CONFIG_FILES=\
	python.xml

PYTHON_FILES=\
	rms_coil_combine.py \
	kspaceandimage.py \
	remove_2x_oversampling.py \
	accumulate_and_recon.py \
	GadgetronXML.py \
	image_viewer.py 

LIBSOURCES=\
	PythonGadget.cpp \
	GadgetReference.cpp \
	GadgetronPythonMRI.cpp \
	PythonCommunicator.cpp 

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

LIBFILE=lib$(DLLNAME).$(DLLEXTENSION)

LIBLDFLAGS += -lfftw3 -lfftw3f -lACE -ltinyxml -lboost_python

PYTHONINC =  $(shell python-config --includes)
PYTHONLIB = $(shell python-config --ldflags)

CXXFLAGS += $(PYTHONINC) -DNPY_NO_DEPRECATED_API
LIBLDFLAGS += $(PYTHONLIB)

all: $(LIBFILE) GadgetronPythonMRI.so

$(LIBFILE): $(LIBOBJECTS) Makefile python_mods
	$(CXX) $(LIBLDFLAGS) $(LIBOBJECTS) -o $@
	mkdir -p $(GADGETRONHOME)/lib
	mkdir -p $(GADGETRONHOME)/config
	cp $(LIBFILE) $(GADGETRONHOME)/lib/
	cp $(CONFIG_FILES) $(GADGETRONHOME)/config/

GadgetronPythonMRI.so: GadgetronPythonMRI.o GadgetReference.o
	$(CXX) $(LIBLDFLAGS) GadgetronPythonMRI.o GadgetReference.o -o GadgetronPythonMRI.so
	cp GadgetronPythonMRI.so  $(GADGETRONHOME)/lib/

.cpp.o:
	$(CXX) $(CXXFLAGS) -DGADGETS_BUILD_DLL $< -o $@

python_mods: $(PYTHON_FILES)
	cp $(PYTHON_FILES) $(GADGETRONHOME)/lib/	

clean:
	rm -rf *~
	rm -rf $(LIBFILE)
	rm -rf *.o
	rm -rf *.cplx
	rm -rf *.real

