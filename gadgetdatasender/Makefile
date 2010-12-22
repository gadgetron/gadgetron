#Gadgetron Makefile

UNAME := $(shell uname)

EXESOURCES=main.cpp siemensraw.cpp
LIBSOURCES=

EXEOBJECTS=$(EXESOURCES:.cpp=.o)

CXX=g++
EXELDFLAGS= -lACE 
CXXFLAGS=-c -Wall -I. #-DACE_NTRACE=0

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
