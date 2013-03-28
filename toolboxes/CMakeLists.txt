add_subdirectory(core)

IF(CUDA_FOUND OR ARMADILLO_FOUND)
  add_subdirectory(operators)
  add_subdirectory(solvers)
ENDIF(CUDA_FOUND OR ARMADILLO_FOUND)

add_subdirectory(mri)

IF(CUDA_FOUND)
  add_subdirectory(nfft)
ENDIF(CUDA_FOUND)

IF (ACE_FOUND)
  add_subdirectory(gadgettools)
ENDIF (ACE_FOUND)
