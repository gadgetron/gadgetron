#ifndef _TESTINGS_H_
#define _TESTINGS_H_

#ifndef min
#define min(a,b)  (((a)<(b))?(a):(b))
#endif
#ifndef max
#define max(a,b)  (((a)<(b))?(b):(a))
#endif

#define TESTING_CUDA_INIT()						\
  CUdevice  dev;							\
  CUcontext context;							\
  if( CUDA_SUCCESS != cuInit( 0 ) ) {					\
    fprintf(stderr, "CUDA: Not initialized\n" ); exit(-1);		\
  }									\
  if( CUDA_SUCCESS != cuDeviceGet( &dev, 0 ) ) {			\
    fprintf(stderr, "CUDA: Cannot get the device\n"); exit(-1);		\
  }									\
  if( CUDA_SUCCESS != cuCtxCreate( &context, 0, dev ) ) {		\
    fprintf(stderr, "CUDA: Cannot create the context\n"); exit(-1);	\
  }									\
  if( CUDA_SUCCESS != cublasInit( ) ) {					\
    fprintf(stderr, "CUBLAS: Not initialized\n"); exit(-1);		\
  }									\
  printout_devices( );


#define TESTING_CUDA_FINALIZE()			\
  cuCtxDetach( context );			\
  cublasShutdown();
     

#define TESTING_MALLOC(ptr, type, size)					\
  ptr = (type*)malloc((size) * sizeof(type));				\
  if (ptr == 0) {							\
    fprintf (stderr, "!!!! Malloc failed for: %s\n", #ptr );		\
    exit(-1);								\
  }

#define TESTING_HOSTALLOC(ptr, type, size)				\
  if( CUBLAS_STATUS_SUCCESS != cudaMallocHost( (void**)&ptr, (size)*sizeof(type) ) ) { \
    fprintf (stderr, "!!!! cudaMallocHost failed for: %s\n", #ptr );	\
    exit(-1);								\
  }

#define TESTING_DEVALLOC(ptr, type, size)				\
  if( CUBLAS_STATUS_SUCCESS != cublasAlloc( size, sizeof(type), (void**)&ptr) ) { \
    fprintf (stderr, "!!!! cublasAlloc failed for: %s\n", #ptr );	\
    exit(-1);								\
  }


#define TESTING_FREE(ptr)			\
  free(ptr);

#define TESTING_HOSTFREE(ptr)			\
  cudaFreeHost( ptr );

#define TESTING_DEVFREE(ptr)			\
  cublasFree( ptr );

#endif /* _TESTINGS_H_ */
