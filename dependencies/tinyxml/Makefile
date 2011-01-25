INSTALL_PREFIX=/usr/local

MACHINE   := $(shell uname -m)
KERNEL    := $(shell uname -s)

DLLEXTENSION=so
ifeq ($(KERNEL), Darwin)
DLLEXTENSION=dylib
endif


HEADERS=\
	tinyxml.h

LIBSOURCES=\
	tinyxml.cpp \
	tinyxmlerror.cpp \
	tinyxmlparser.cpp

LIBOBJECTS=$(LIBSOURCES:.cpp=.o)

CXX=g++
ifeq ($(KERNEL), Darwin)
CXX=g++ -m32 -arch i386
endif

LIBLDFLAGS=-shared

CXXFLAGS=-c -fPIC -Wall -I. 


LIBFILE=libtinyxml.$(DLLEXTENSION)

all: $(LIBFILE) 

$(LIBFILE): $(LIBOBJECTS) Makefile
	$(CXX) $(LIBLDFLAGS) $(LIBOBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -rf *~
	rm -rf $(LIBFILE)
	rm -rf *.o

install:
	cp $(HEADERS) $(INSTALL_PREFIX)/include/
	cp $(LIBFILE) $(INSTALL_PREFIX)/lib/
