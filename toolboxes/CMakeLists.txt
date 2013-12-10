add_subdirectory(core)

IF(CUDA_FOUND OR ARMADILLO_FOUND)
  add_subdirectory(operators)
  add_subdirectory(solvers)
ENDIF(CUDA_FOUND OR ARMADILLO_FOUND)

add_subdirectory(mri)
add_subdirectory(nfft)
add_subdirectory(registration)

IF (ACE_FOUND AND XSD_FOUND)
  add_subdirectory(gadgettools)
ENDIF (ACE_FOUND AND XSD_FOUND)

IF (MKL_FOUND)
    add_subdirectory(gtplus)
ENDIF (MKL_FOUND)