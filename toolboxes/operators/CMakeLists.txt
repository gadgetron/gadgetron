

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
  DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)

add_library(pingvin_toolbox_operator INTERFACE)
target_include_directories(pingvin_toolbox_operator
        INTERFACE
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

install(TARGETS pingvin_toolbox_operator
        EXPORT pingvin-export)

add_subdirectory(cpu)