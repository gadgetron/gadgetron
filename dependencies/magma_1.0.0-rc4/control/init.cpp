/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010
*/

#include <stdio.h>
#include <string.h>
#include <pthread.h>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>

#include <quark.h>

#include "magma.h"

typedef struct {
  int tid;
  void *params;
} t_params;

extern "C" magma_context *
magma_init( void *params, void* (*func)(void *a), magma_int_t nthread, magma_int_t ncpu, magma_int_t ngpu,
	    magma_int_t argc, char **argv)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======
    This function initializes the hardware context to be used for
    subsequent calls to routines in the MAGMA library.

    Arguments
    =========
    NCPU   (input) INTEGER
            Number of CPU cores to be used in the computations.

    NGPU    (input) INTEGER
            Number of GPU cores to be used in the computations.
    ===================================================================== */

    t_params **tp = (t_params**)malloc(sizeof(t_params*)*nthread);

    pthread_t *thread;

    magma_int_t i;
    magma_context *context;
    context  = (magma_context *)malloc(sizeof(magma_context));

    if (nthread > 0) {
      thread = (pthread_t*)malloc(sizeof(pthread_t)*nthread);

      for (i = 0; i < nthread; i++){
        tp[i] = (t_params*)malloc(sizeof(t_params));
		tp[i]->params = params;
		tp[i]->tid = i;
        pthread_create(&thread[i], NULL, func, (void *)tp[i]);
      }
    }


    if (ncpu <= 1)
      ncpu = 1;

    if (ngpu <= 0)
      ngpu = 0;

    context->num_cores = ncpu;
    context->num_gpus  = ngpu; 

    if (ncpu > 1)
      {
	    /* Initialize the QUARK scheduler */
        context->quark = QUARK_New(ncpu);
      }

    if (ngpu > 1)
      {
        printf("The requested number of GPUs is not yet supported.\n\n");
        printf("The number of GPUs set to one.\n\n");
        context->num_gpus  = 1;
      }
    
    if (ngpu == 1)
      {
	CUdevice  dev;
	context->gpu_context = (CUcontext *)malloc(ngpu * sizeof(CUcontext));
	
	/* For now we use by default device 0, always */
	if( CUDA_SUCCESS != cuInit( 0 ) ) {
	  fprintf(stderr, "CUDA: Not initialized\n" ); 
	  exit(-1);
	}
	if( CUDA_SUCCESS != cuDeviceGet( &dev, 0 ) ) {
	  fprintf(stderr, "CUDA: Cannot get the device\n");
	  exit(-1);
	} 
	if( CUDA_SUCCESS != cuCtxCreate( &context->gpu_context[0], 0, dev ) ) {
	  fprintf(stderr, "CUDA: Cannot create the context\n");
	  exit(-1);
	}
	if( CUDA_SUCCESS != cublasInit( ) ) {
	  fprintf(stderr, "CUBLAS: Not initialized\n");
	  exit(-1);
	}
	printout_devices( );
      }
    
    context->nb = -1;
    for(i = 1; i<argc; i++)
      if (strcmp("-b", argv[i])==0)
	context->nb = atoi(argv[++i]);
    
    return context;
}


extern "C" void
magma_finalize( magma_context *cntxt)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======
    This function finalizes the MAGMA hardware context.

    Arguments
    =========
    CNTXT  (input) MAGMA_CONTEXT
           Pointer to the MAGMA hardware context to be closed
    ===================================================================== */

  if (cntxt->num_cores > 1)
    /* Shut down the QUARK scheduler */
    QUARK_Delete(cntxt->quark);

  if (cntxt->num_gpus == 1)
    {
      /* Shutdown CUDA and CUBLAS*/
      cuCtxDetach( cntxt->gpu_context[0] );
      cublasShutdown();

      free(cntxt->gpu_context);
    }

  free(cntxt);
}
