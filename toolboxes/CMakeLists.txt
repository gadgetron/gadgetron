add_subdirectory(ndarray)
add_subdirectory(gadgettools)
add_subdirectory(hostutils)

if (BLAS_FOUND)
 	add_subdirectory(linalg)
endif(BLAS_FOUND)

if (CUDA_FOUND)
   add_subdirectory(gpucore)
   add_subdirectory(gpunfft)
   add_subdirectory(gpupmri)
   add_subdirectory(solvers)
endif(CUDA_FOUND)


