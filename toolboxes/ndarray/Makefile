GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc

LIBHEADERS=\
	NDArray.h \
	hoNDArray.h \
	hoNDArray_fileio.h \
	FFT.h \
	gadgetron_export.h

install: $(LIBHEADERS)
	mkdir -p $(GADGETRONHOME)/inc
	cp $(LIBHEADERS) $(GADGETRONHOME)/inc/

clean:
	rm -rf *~
	rm -rf $(LIBFILE)
	rm -rf *.o
	rm -rf *.cplx
	rm -rf *.real

