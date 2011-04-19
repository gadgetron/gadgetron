GADGETRONHOME=.
include $(GADGETRONHOME)/Makefile.inc

DIRS=\
	toolboxes/gadgettools \
	toolboxes/ndarray \
	toolboxes/gpucore \
	toolboxes/gpuNFFT \
	toolboxes/gpuCSM \
	toolboxes/gpucg \
	apps/gputest \
	apps/gadgetron \
	apps/gadgetdatasender \
	apps/htgrappa \
	gadgets/core \
	gadgets/gpucg

all: $(DIRS)

$(DIRS): force_look
	cd $@ ; $(MAKE) $(MFLAGS)

clean :
	rm -rf inc
	rm -rf lib
	rm -rf bin
	rm -rf *~
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look :
	true