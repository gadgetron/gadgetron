#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 1.0) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      November 2010
#//////////////////////////////////////////////////////////////////////////////

MAGMA_DIR = .
include ./Makefile.internal

all: lib test

lib: libmagma libmagmablas

clean: cleanall

libmagma:
	( cd control     && $(MAKE) )
	( cd src         && $(MAKE) )

libmagmablas:
	( cd magmablas   && $(MAKE) )

libquark:
	( cd quark       && $(MAKE) )

test:
	( cd testing/lin && $(MAKE) )
	( cd testing     && $(MAKE) )

clean:
	( cd control     && $(MAKE) clean )
	( cd src         && $(MAKE) clean )
	( cd testing     && $(MAKE) clean )
	( cd testing/lin && $(MAKE) clean )
	( cd magmablas   && $(MAKE) clean ) 
#	( cd quark       && $(MAKE) clean )

cleanall:
	( cd control     && $(MAKE) cleanall )
	( cd src         && $(MAKE) cleanall )
	( cd testing     && $(MAKE) cleanall )
	( cd testing/lin && $(MAKE) cleanall )
	( cd magmablas   && $(MAKE) cleanall ) 
	( cd lib && rm -f *.a )
#	( cd quark       && $(MAKE) cleanall )


dir:
	mkdir -p $(prefix)
	mkdir -p $(prefix)/include
	mkdir -p $(prefix)/lib
	mkdir -p $(prefix)/lib/pkgconfig

install: lib dir
#       MAGMA
	cp $(MAGMA_DIR)/include/*.h  $(prefix)/include
	cp $(LIBMAGMA)               $(prefix)/lib
	cp $(LIBMAGMABLAS)           $(prefix)/lib
#       QUARK
#	cp $(QUARKDIR)/include/quark.h             $(prefix)/include
#	cp $(QUARKDIR)/include/quark_unpack_args.h $(prefix)/include
#	cp $(QUARKDIR)/include/icl_hash.h          $(prefix)/include
#	cp $(QUARKDIR)/include/icl_list.h          $(prefix)/include
#	cp $(QUARKDIR)/lib/libquark.a              $(prefix)/lib
#       pkgconfig
	cat $(MAGMA_DIR)/lib/pkgconfig/magma.pc | \
	    sed -e s+\__PREFIX+"$(prefix)"+     | \
	    sed -e s+\__LIBEXT+"$(LIBEXT)"+       \
	    > $(prefix)/lib/pkgconfig/magma.pc

