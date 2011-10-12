GADGETRONHOME=.
include $(GADGETRONHOME)/Makefile.inc

DIRS=\
	toolboxes/gadgettools \
	toolboxes/ndarray \
	toolboxes/hostutils \
	toolboxes/gpucore \
	toolboxes/gpunfft \
	toolboxes/solvers \
	toolboxes/gpupmri \
	toolboxes/gpuct \
	apps/standalone/gpu/gputest \
	apps/standalone/gpu/MRI/nfft/2d \
	apps/standalone/gpu/MRI/nfft/ms2d \
	apps/standalone/gpu/MRI/sense/noncartesian/radial/2d_golden_ratio \
	apps/standalone/gpu/MRI/sense/noncartesian/radial/2d_golden_ratio_kt \
	apps/standalone/gpu/CT/parallel_beam/2d \
	apps/standalone/gpu/CT/parallel_beam/2d+t \
	apps/gadgetron \
	apps/gadgetdatasender \
	gadgets/core \
	gadgets/grappa \
	gadgets/python \
	gadgets/cgsense
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