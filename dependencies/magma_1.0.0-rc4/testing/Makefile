#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 1.0) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      November 2010
#//////////////////////////////////////////////////////////////////////////////

MAGMA_DIR = ..
include ../Makefile.internal

ZSRC  = testing_zgemm.cpp	\
	testing_zhemv.cpp	\
	testing_zgemv.cpp	\
	testing_zsymv.cpp	\
	\
	testing_zpotrf.cpp	\
	testing_zgetrf.cpp 	\
	testing_zgeqrf.cpp 	\
	testing_zgeqlf.cpp 	\
	testing_zgelqf.cpp 	\
	\
	testing_zpotrf_gpu.cpp 	\
	testing_zgetrf_gpu.cpp 	\
	testing_zgesv_gpu.cpp	\
	testing_zgeqrf_gpu.cpp 	\
        testing_zgeqrs_gpu.cpp 	\
	testing_zungqr_gpu.cpp	\
	\
	testing_zgebrd.cpp	\
        testing_zgehrd.cpp 	\
        testing_zhetrd.cpp 	\
        testing_zgeev.cpp       \
        testing_zgesvd.cpp      \
        testing_zheevd.cpp      \
        testing_dsyevd.cpp	\
	\
	testing_zcposv_gpu.cpp	\
	testing_zcgesv_gpu.cpp	\
	testing_zcgeqrsv_gpu.cpp

#testing_zswap.cpp

ifeq ($(FORT), pgfortran)
	FOBJ := fortran_thunking.o
	ZSRCF = testing_zgetrf_gpu_f.cuf
else
	FOBJ := fortran.o
	ZSRCF = testing_zgetrf_gpu_f.f
endif

ZSRC += $(ZSRCF) 
CSRC = testing_cgemm.cpp testing_chemv.cpp testing_cgemv.cpp testing_csymv.cpp testing_cpotrf.cpp testing_cgetrf.cpp testing_cgeqrf.cpp testing_cgeqlf.cpp testing_cgelqf.cpp testing_cpotrf_gpu.cpp testing_cgetrf_gpu.cpp testing_cgesv_gpu.cpp testing_cgeqrf_gpu.cpp testing_cgeqrs_gpu.cpp testing_cungqr_gpu.cpp testing_cgebrd.cpp testing_cgehrd.cpp testing_chetrd.cpp testing_cgeev.cpp testing_cgesvd.cpp testing_cheevd.cpp testing_cgetrf_gpu_f.f
DSRC = testing_dgemm.cpp testing_dsymv.cpp testing_dgemv.cpp testing_dpotrf.cpp testing_dgetrf.cpp testing_dgeqrf.cpp testing_dgeqlf.cpp testing_dgelqf.cpp testing_dpotrf_gpu.cpp testing_dgetrf_gpu.cpp testing_dgesv_gpu.cpp testing_dgeqrf_gpu.cpp testing_dgeqrs_gpu.cpp testing_dorgqr_gpu.cpp testing_dgebrd.cpp testing_dgehrd.cpp testing_dsytrd.cpp testing_dgeev.cpp testing_dgesvd.cpp testing_dsposv_gpu.cpp testing_dsgesv_gpu.cpp testing_dsgeqrsv_gpu.cpp testing_dgetrf_gpu_f.f
SSRC = testing_sgemm.cpp testing_ssymv.cpp testing_sgemv.cpp testing_spotrf.cpp testing_sgetrf.cpp testing_sgeqrf.cpp testing_sgeqlf.cpp testing_sgelqf.cpp testing_spotrf_gpu.cpp testing_sgetrf_gpu.cpp testing_sgesv_gpu.cpp testing_sgeqrf_gpu.cpp testing_sgeqrs_gpu.cpp testing_sorgqr_gpu.cpp testing_sgebrd.cpp testing_sgehrd.cpp testing_ssytrd.cpp testing_sgeev.cpp testing_sgesvd.cpp testing_ssyev.cpp testing_sgetrf_gpu_f.f

ALLSRC  = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJF  = $(ALLSRC:.f=.o)
ALLOBJF := $(ALLOBJF:.cuf=.o)
ALLOBJ   = $(ALLOBJF:.cpp=.o)
EXE      = $(ALLOBJ:.o=)

all: $(EXE)

$(EXE): $(LIBMAGMA) $(LIBMAGMABLAS)

clean:
	rm -f *.o *~

cleanall: clean
	rm -f $(EXE)

ifeq ($(FORT), pgfortran)
testing_%.o : testing_%.cuf
	$(FORT) $(OPTS) $(INC) -c $< -o $@ 
else
testing_%.o : testing_%.f
	$(FORT) $(OPTS) $(INC) -c $< -o $@ 
endif

%.o: %.cpp
	$(CC) $(OPTS) $(INC) -c $< -o $@

#
# MAGMA Fortran module. Once finalized, the generated magma.mod
# has to be moved in the include directory
#
testing_%_f : testing_%_f.o $(FOBJ)
	$(FORT) $(OPTS) $(LDOPTS) $< $(FOBJ) -o $@ lin/liblapacktest.a -L../lib \
          -lcuda -lmagma -lmagmablas -lmagma $(LIBDIR) $(LIB) 

testing_% : testing_%.o
	$(CC) $(OPTS) $(LDOPTS) $(NOMAIN) $< $(OBJ) -o $@ lin/liblapacktest.a -L../lib \
          -lcuda -lmagma -lmagmablas -lmagma $(LIBDIR) $(LIB) 

#
# These two files are provided by cuda
#
fortran.o : $(CUDADIR)/src/fortran.c
	$(CC) $(OPTS) $(INC) -c $< -o $@

fortran_thunking.o : $(CUDADIR)/src/fortran_thunking.c
	$(CC) $(OPTS) $(INC) -c $< -o $@

