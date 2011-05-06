#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 1.0) --
#      Univ. of Tennessee
#      Univ. of California Berkeley
#      November 2010
#//////////////////////////////////////////////////////////////////////////////

MAGMA_DIR = ..
include ../Makefile.internal

INC += -I$(MAGMA_DIR)/control

ZHDR  = commonblas_z.h
CHDR = commonblas_c.h
DHDR = commonblas_d.h
SHDR = commonblas_s.h

ALLHDR= $(ZHDR) $(CHDR) $(DHDR) $(SHDR) \
        commonblas.h

ZSRC  = zauxiliary.cu	\
	zpermute.cu	\
	zpermute-v2.cu	\
	zswap.cu	\
	zswapblk.cu	\
	zlacpy.cu	\
	zlange.cu       \
	zlanhe.cu	\
	ztranspose.cu 	\
	ztranspose-v2.cu\
	zinplace_transpose.cu	\
	\
	zclaswp.cu 	\
        zcaxpycp.cu     \
	zlag2c.cu       \
	clag2z.cu       \
	zlat2c.cu       \
	zsymv_fermi.cu	\
	zsymv_tesla.cu	\
	csymv_tesla.cu	\
	zhemv_fermi.cu	\
	zhemv_tesla.cu	\
	chemv_tesla.cu

CSRC = cauxiliary.cu cpermute.cu cpermute-v2.cu cswap.cu cswapblk.cu clacpy.cu clange.cu clanhe.cu ctranspose.cu ctranspose-v2.cu cinplace_transpose.cu csymv_fermi.cu chemv_fermi.cu
DSRC = dauxiliary.cu dpermute.cu dpermute-v2.cu dswap.cu dswapblk.cu dlacpy.cu dlange.cu dlansy.cu dtranspose.cu dtranspose-v2.cu dinplace_transpose.cu dslaswp.cu dsaxpycp.cu dlag2s.cu slag2d.cu dlat2s.cu dsymv_fermi.cu dsymv_tesla.cu
SSRC = sauxiliary.cu spermute.cu spermute-v2.cu sswap.cu sswapblk.cu slacpy.cu slange.cu slansy.cu stranspose.cu stranspose-v2.cu sinplace_transpose.cu ssymv_fermi.cu ssymv_tesla.cu

FERMISRC =	sgemv_fermi.cu  \
                dgemv_fermi.cu  \
                sgemm_fermi.cu 	\
                cgemv_fermi.cu 	\
                zgemv_fermi.cu 	\
                sgemm_fermi80.cu\
                sgemm_fermi64.cu\
	        dgemm_fermi.cu 	\
                \
                strsm_tesla.cu 	\
                dtrsm_tesla.cu	\
                \
                ssyr2k.cu       \

TESLASRC =      sgemv_tesla.cu	\
	        dgemv_tesla.cu 	\
	        cgemv_tesla.cu 	\
	        zgemv_tesla.cu 	\
	        gemv32_tesla.cu \
                strsm_tesla.cu  \
                dtrsm_tesla.cu  \
                \
	        sgemm_tesla.cu 					\
		sgemm_tesla_a_0.cu 				\
		sgemm_tesla_ab_0.cu 				\
		sgemm_tesla_N_N_64_16_16_16_4.cu 		\
		sgemm_tesla_N_N_64_16_16_16_4_special.cu 	\
		sgemm_tesla_N_T_64_16_4_16_4.cu 		\
		sgemm_tesla_T_N_32_32_8_8_8.cu 			\
		sgemm_tesla_T_T_64_16_16_16_4.cu 		\
		\
		dgemm_tesla.cu 					\
		dgemm_tesla_a_0.cu 				\
		dgemm_tesla_ab_0.cu 				\
		dgemm_tesla_N_N_64_16_16_16_4.cu 		\
		dgemm_tesla_N_N_64_16_16_16_4_special.cu 	\
		dgemm_tesla_N_T_64_16_4_16_4.cu 		\
		dgemm_tesla_T_N_32_32_8_8_8.cu 			\
		dgemm_tesla_T_T_64_16_16_16_4.cu 		\
		dgemm_tesla_T_T_64_16_16_16_4_v2.cu             \
                \
                ssyr2k.cu                                       \

# 		zhemv_tesla.cu

ALLSRC= $(ZSRC) $(CSRC) $(DSRC) $(SSRC) \
        dgemv_MLU.cu

ifeq (${GPU_TARGET}, 0)
ALLSRC += $(TESLASRC)
GPUOPT = ${TESLAOPT}
endif

ifeq (${GPU_TARGET}, 1)
ALLSRC += $(FERMISRC)
GPUOPT = ${FERMIOPT}
endif

ALLOBJ = $(ALLSRC:.cu=.cu_o) 

all: $(LIBMAGMABLAS) 

$(LIBMAGMABLAS): $(ALLHDR) $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

%.o : %.cpp 
	$(CC) $(OPTS) $(INC) -c $< -o $@

%.cu_o: %.cu
	$(NVCC) $(NVOPTS) $(GPUOPT) $(INC) -c $< -o $@

clean:
	rm -f *.cu_o *.o *~ *.linkinfo

cleanall: clean
	rm -f *.a
	rm -f $(LIBMAGMABLAS)
