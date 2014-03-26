add_subdirectory(core)

add_subdirectory(operators)
add_subdirectory(solvers)

add_subdirectory(fft)
add_subdirectory(nfft)

add_subdirectory(mri)
add_subdirectory(ct)

add_subdirectory(registration)

IF (ACE_FOUND AND XSD_FOUND AND XERCESC_FOUND)
  add_subdirectory(gadgettools)
ENDIF()

IF (MKL_FOUND AND FFTW3_FOUND AND ISMRMRD_FOUND)
  add_subdirectory(gtplus)
ENDIF()
