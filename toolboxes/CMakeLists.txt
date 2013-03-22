add_subdirectory(core)

IF (ACE_FOUND)
  add_subdirectory(gadgettools)
ENDIF (ACE_FOUND)

if (CUDA_FOUND)
#   add_subdirectory(gpucore)
#   add_subdirectory(gpunfft)
#   add_subdirectory(gpupmri)
#   add_subdirectory(solvers)
#   add_subdirectory(operators)
endif(CUDA_FOUND)
