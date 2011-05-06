/**
 *
 * @file context.h
 *
 *  MAGMA (version 1.0) --
 *  Univ. of Tennessee, Knoxville
 *  Univ. of California, Berkeley
 *  Univ. of Colorado, Denver
 *  November 2010
 *
 **/

#ifndef _MAGMA_CONTEXT_H_
#define _MAGMA_CONTEXT_H_


/* ------------------------------------------------------------
 * MAGMA Context
 * --------------------------------------------------------- */
typedef struct magma_context_s
{
  /* Number of CPU core in this context */
  magma_int_t num_cores;

  /* Number of GPUs in this context */
  magma_int_t num_gpus;

  /* GPU contexts */
  CUcontext *gpu_context;

  /* QUARK scheduler */
  Quark *quark;

  /* Block size, internally used for some algorithms */
  magma_int_t nb;

  /* Pointer to other global algorithm-dependent parameters */ 
  void *params;

} magma_context_t;



magma_context_t * magma_init(void *, void* (*func)(void *a), magma_int_t nthread, magma_int_t ncpu, magma_int_t ngpu,
                magma_int_t argc, char **argv);
void magma_finalize(magma_context_t * cntxt);

void auto_tune(char algorithm, char precision, magma_int_t ncores, magma_int_t ncorespsocket,
               magma_int_t m, magma_int_t n, magma_int_t *nb, magma_int_t *ob, magma_int_t *ib,
               magma_int_t *nthreads, magma_int_t *nquarkthreads);  


#endif /* _MAGMA_CONTEXT_H_ */
