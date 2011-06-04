GADGETRONHOME=.
include $(GADGETRONHOME)/Makefile.inc

DIRS=\
	toolboxes/gadgettools \
	toolboxes/ndarray \
	toolboxes/hostutils \
	toolboxes/gpucore \
	toolboxes/gpunfft \
	toolboxes/gpucg \
	toolboxes/gpupmri \
	apps/MRI/gputest \
	apps/MRI/standalone/radial_sense \
	apps/gadgetron \
	apps/gadgetdatasender \
	gadgets/core \
	gadgets/grappa
#	gadgets/gpucg \

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