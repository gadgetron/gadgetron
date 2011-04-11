GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

HEADERS=\
	GadgetContainerMessage.h \
	Gadget.h \
	GadgetServerAcceptor.h \
	GadgetStreamController.h \
	Gadgetron.h \
	GadgetronExport.h \
	GadgetMessageInterface.h \

EXESOURCES=\
	main.cpp \
	GadgetStreamController.cpp \
	GadgetServerAcceptor.cpp

EXEOBJECTS=$(EXESOURCES:.cpp=.o)
LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

EXECUTABLE=gadgetron

all: $(EXECUTABLE)

$(EXECUTABLE): $(EXEOBJECTS) Makefile
	$(CXX) $(APPLDFLAGS) $(EXEOBJECTS) -o $@
	mkdir -p $(GADGETRONHOME)/inc
	mkdir -p $(GADGETRONHOME)/bin
	cp $(HEADERS) $(GADGETRONHOME)/inc/
	cp $(EXECUTABLE) $(GADGETRONHOME)/bin/

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -rf *~
	rm -rf $(EXECUTABLE)
	rm -rf *.o
	rm -rf *.cplx
	rm -rf *.real