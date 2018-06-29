include_directories(
  ${Boost_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_BINARY_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  )

set(operator_header_files 
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/generalOperator.h
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/linearOperator.h
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/identityOperator.h
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/diagonalOperator.h
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/diagonalSumOperator.h
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/encodingOperatorContainer.h
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/multiplicationOperatorContainer.h
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/FFTOperator.h
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/imageOperator.h
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/encodedImageOperator.h
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/partialDerivativeOperator.h
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/convolutionOperator.h
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/laplaceOperator.h
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/downsampleOperator.h
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/upsampleOperator.h
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/tvPicsOperator.h
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/subsetOperator.h
    )
  
install(FILES ${operator_header_files}
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

if (ARMADILLO_FOUND)
  add_subdirectory(cpu)
endif ()

if (CUDA_FOUND)
  add_subdirectory(gpu)
endif ()
