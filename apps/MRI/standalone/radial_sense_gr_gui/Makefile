GADGETRONHOME=../../../..
include $(GADGETRONHOME)/Makefile.inc

EXECUTABLE = radial_sense_gr_gui

UI_ELEMENTS = reconBaseWidget.ui radialSenseAppBaseMainWidget.ui
UI_SOURCES = radialSenseAppMainWidget.cpp reconWidget.cpp GLReconWidget.cpp
MOC_HEADERS = radialSenseAppMainWidget.h reconWidget.h GLReconWidget.h 
MOC_SOURCES = moc_radialSenseAppMainWidget.cpp moc_reconWidget.cpp moc_GLReconWidget.cpp

HEADERFILES =

APPSOURCES = main.cpp $(UI_SOURCES)

CUDASOURCES =

UIOBJECTS = $(UI_ELEMENTS:.ui=.h)
MOCOBJECTS = $(MOC_HEADERS:.h=.cpp)
CUDAOBJECTS = $(CUDASOURCES:.cu=.cuo)
APPOBJECTS = $(APPSOURCES:.cpp=.o)
MOCSOBJECTS = $(MOC_SOURCES:.cpp=.o)

# cuda sdk root
CUDA_SDK_DIR ?= $(HOME)/NVIDIA_GPU_Computing_SDK

# automagically find cuda libs on osx (without setting DYLD_LIBRARY_PATH)
ifeq ($(KERNEL), Darwin)
CUDA_INSTALL_PATH ?= /usr/local/cuda
APPLDFLAGS += -Xlinker -rpath $(CUDA_INSTALL_PATH)/lib
endif

# Qt libs
ifeq ($(KERNEL), Darwin)
APPLDFLAGS += -framework QtGui -framework QtCore -framework QtOpenGL -framework OpenGL -lglew
else
CXXFLAGS += -I/usr/include/qt4
APPLDFLAGS += -L$(CUDA_SDK_DIR)/C/common/lib/linux -lGLEW_x86_64 -lQtGui -lQtCore -lQtOpenGL
endif

# gadgetron libs
APPLDFLAGS += -lgpucore -lgpunfft -lgpucg -lgpuparallelmri

%.h: %.ui
	uic $< -o ui_$@

%.cpp: %.h
	moc $< -o moc_$@

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

%.cuo: %.cu
	$(NVCC) $(NVCCFLAGS) -o $@ $<

gui: $(UIOBJECTS) $(MOCOBJECTS) $(CUDAOBJECTS) $(APPOBJECTS) $(MOCSOBJECTS)
	$(CXX) $(APPOBJECTS) $(CUDAOBJECTS) $(MOCSOBJECTS) $(CUDALIBRARIES) $(APPLDFLAGS) -o $(EXECUTABLE) 
	mkdir -p $(GADGETRONHOME)/bin
	cp $(EXECUTABLE) $(GADGETRONHOME)/bin/

clean:
	rm -rf *~
	rm -rf $(EXECUTABLE)
	rm -rf $(APPOBJECTS)
	rm -rf $(CUDAOBJECTS)
	rm -rf moc_*
	rm -rf ui_*.h
	touch *.h
