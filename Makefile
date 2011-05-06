GADGETRONHOME=.
include $(GADGETRONHOME)/Makefile.inc

DIRS=\
	toolboxes/gadgettools \
	toolboxes/hostutils \
	toolboxes/ndarray \
	toolboxes/gpucore \
	toolboxes/gpuNFFT \
	toolboxes/gpucg \
	toolboxes/gpuParallelMRI \
	apps/gputest \
	apps/standalone/radial_sense \
	apps/gadgetron \
	apps/gadgetdatasender \
	gadgets/core \
	gadgets/gpucg \
	gadgets/grappa

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