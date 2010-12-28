#Gadgetron Makefile

UNAME := $(shell uname)

EXESOURCES=main.cpp siemensraw.cpp ../gadgettools/ConfigParser.cpp
LIBSOURCES=

EXEOBJECTS=$(EXESOURCES:.cpp=.o)

CXX=g++
EXELDFLAGS= -lACE 
CXXFLAGS=-c -Wall -I. -I../gadgettools #-DACE_NTRACE=0

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
