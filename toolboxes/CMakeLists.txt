add_subdirectory(cpucore)
add_subdirectory(hostutils)

IF (ACE_FOUND)
        add_subdirectory(gadgettools)
ENDIF (ACE_FOUND)

if (BLAS_FOUND)
 	add_subdirectory(linalg)
endif(BLAS_FOUND)

if (CUDA_FOUND)
   add_subdirectory(gpucore)
   add_subdirectory(gpunfft)
   add_subdirectory(gpupmri)
   add_subdirectory(solvers)
endif(CUDA_FOUND)


