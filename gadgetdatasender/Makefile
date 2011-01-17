#Gadgetron Makefile

GADGETRONHOME=..

MACHINE   := $(shell uname -m)
KERNEL    := $(shell uname -s)

EXESOURCES=main.cpp siemensraw.cpp
EXEOBJECTS=$(EXESOURCES:.cpp=.o)

CXX=g++
ifeq ($(KERNEL), Darwin)
CXX=g++ -m32 -arch i386
endif

EXELDFLAGS= -L$(GADGETRONHOME)/lib -lACE -lgadgetrontools 
CXXFLAGS=-c -Wall -I. -I$(GADGETRONHOME)/inc #-DACE_NTRACE=0

EXECUTABLE=gadgetrondatasender

all: $(EXECUTABLE)

$(EXECUTABLE): $(EXEOBJECTS) Makefile
	$(CXX) $(EXELDFLAGS) $(EXEOBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -rf *~
	rm -rf $(EXECUTABLE)
	rm -rf *.o
	rm -rf *.cplx
	rm -rf *.real
