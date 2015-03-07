include_directories(
  ${Boost_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  )

install(FILES 	
  generalOperator.h
  linearOperator.h
  identityOperator.h
  diagonalOperator.h
  diagonalSumOperator.h
  encodingOperatorContainer.h
  multiplicationOperatorContainer.h
  FFTOperator.h
  imageOperator.h
  encodedImageOperator.h
  partialDerivativeOperator.h
  convolutionOperator.h
  laplaceOperator.h
  downsampleOperator.h
  upsampleOperator.h
  tvPicsOperator.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

IF (ARMADILLO_FOUND)
  add_subdirectory(cpu)
ENDIF (ARMADILLO_FOUND)

IF (CUDA_FOUND)
  add_subdirectory(gpu)
ENDIF (CUDA_FOUND)
