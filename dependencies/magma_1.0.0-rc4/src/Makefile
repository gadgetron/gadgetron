#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 1.0) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      November 2010
#//////////////////////////////////////////////////////////////////////////////

MAGMA_DIR = ..
include ../Makefile.internal

INC += -I$(MAGMA_DIR)/control

ZSRC  = zpotrf.cpp	\
        zpotrf_gpu.cpp	\
        zpotrs_gpu.cpp	\
        zposv_gpu.cpp	\
	\
        zgetrf.cpp	\
        zgetrf_gpu.cpp 	\
        zgetrs_gpu.cpp 	\
        zgesv_gpu.cpp 	\
	\
	zgeqrf.cpp	\
	zgeqlf.cpp	\
	zgelqf.cpp	\
        zungqr.cpp	\
	zgeqrf_gpu.cpp  \
	zgeqrf2_gpu.cpp	\
	zgeqrs_gpu.cpp	\
	zgels_gpu.cpp	\
	zlarfb_gpu.cpp	\
        zunmqr_gpu.cpp	\
        zungqr_gpu.cpp	\
	\
	zgebrd.cpp 	\
	zlabrd_gpu.cpp	\
	zgehrd.cpp 	\
        zgehrd2.cpp     \
        zlahru.cpp 	\
        zlahr2.cpp 	\
        zunghr.cpp      \
        zhetrd.cpp 	\
        zlatrd.cpp 	\
        zunmtr.cpp 	\
        zunmqr.cpp 	\
        \
        zgeev.cpp       \
        dgeev.cpp       \
        zheevd.cpp      \
        dsyevd.cpp      \
        zgesvd.cpp	\
	dgesvd.cpp	\
	\
	zcposv_gpu.cpp	\
	zcgesv_gpu.cpp	\
	zcgetrs_gpu.cpp	\
	zcgeqrsv_gpu.cpp\
	\
	zgessm_gpu.cpp	\
	zssssm_gpu.cpp
#	zgetfl_gpu.cpp	\
#	ztstrf_gpu.cpp	\

CSRC = cpotrf.cpp cpotrf_gpu.cpp cpotrs_gpu.cpp cposv_gpu.cpp cgetrf.cpp cgetrf_gpu.cpp cgetrs_gpu.cpp cgesv_gpu.cpp cgeqrf.cpp cgeqlf.cpp cgelqf.cpp cungqr.cpp cgeqrf_gpu.cpp cgeqrf2_gpu.cpp cgeqrs_gpu.cpp cgels_gpu.cpp clarfb_gpu.cpp cunmqr_gpu.cpp cungqr_gpu.cpp cgebrd.cpp clabrd_gpu.cpp cgehrd.cpp cgehrd2.cpp clahru.cpp clahr2.cpp cunghr.cpp chetrd.cpp clatrd.cpp cunmtr.cpp cunmqr.cpp cgeev.cpp cheevd.cpp cgesvd.cpp cgessm_gpu.cpp cssssm_gpu.cpp
DSRC = dpotrf.cpp dpotrf_gpu.cpp dpotrs_gpu.cpp dposv_gpu.cpp dgetrf.cpp dgetrf_gpu.cpp dgetrs_gpu.cpp dgesv_gpu.cpp dgeqrf.cpp dgeqlf.cpp dgelqf.cpp dorgqr.cpp dgeqrf_gpu.cpp dgeqrf2_gpu.cpp dgeqrs_gpu.cpp dgels_gpu.cpp dlarfb_gpu.cpp dormqr_gpu.cpp dorgqr_gpu.cpp dgebrd.cpp dlabrd_gpu.cpp dgehrd.cpp dgehrd2.cpp dlahru.cpp dlahr2.cpp dorghr.cpp dsytrd.cpp dlatrd.cpp dormtr.cpp dormqr.cpp dsposv_gpu.cpp dsgesv_gpu.cpp dsgetrs_gpu.cpp dsgeqrsv_gpu.cpp dgessm_gpu.cpp dssssm_gpu.cpp
SSRC = spotrf.cpp spotrf_gpu.cpp spotrs_gpu.cpp sposv_gpu.cpp sgetrf.cpp sgetrf_gpu.cpp sgetrs_gpu.cpp sgesv_gpu.cpp sgeqrf.cpp sgeqlf.cpp sgelqf.cpp sorgqr.cpp sgeqrf_gpu.cpp sgeqrf2_gpu.cpp sgeqrs_gpu.cpp sgels_gpu.cpp slarfb_gpu.cpp sormqr_gpu.cpp sorgqr_gpu.cpp sgebrd.cpp slabrd_gpu.cpp sgehrd.cpp sgehrd2.cpp slahru.cpp slahr2.cpp sorghr.cpp ssytrd.cpp slatrd.cpp sormtr.cpp sormqr.cpp sgeev.cpp ssyev.cpp sgesvd.cpp sgessm_gpu.cpp sssssm_gpu.cpp

ALLSRC= $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ  = $(ALLSRC:.cpp=.o)

all: $(LIBMAGMA)

$(LIBMAGMA): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

.cpp.o:
	$(CC) $(OPTS) $(INC) -c $< -o $@

clean:
	rm -f *.o *~ *.linkinfo

cleanall: clean
	rm -f *.a
	rm -f $(LIBMAGMA)
