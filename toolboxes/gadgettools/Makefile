GADGETRONHOME=../..
include $(GADGETRONHOME)/Makefile.inc


LIBHEADERS=\
	GadgetSocketReceiver.h \
	GadgetSocketSender.h \
	GadgetXml.h

install: $(LIBHEADERS)
	mkdir -p $(GADGETRONHOME)/inc
	cp $(LIBHEADERS) $(GADGETRONHOME)/inc/

clean:
	rm -rf *~
	rm -rf *.o
	rm -rf *.cplx
	rm -rf *.real

