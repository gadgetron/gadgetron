if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_CPUOPERATORS__)
endif ()

include_directories(
    ${ISMRMRD_INCLUDE_DIR}
    ${ARMADILLO_INCLUDE_DIRS}
    ${FFTW3_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/operators
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/dwt/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/application
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/dissimilarity
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/register
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/solver
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/transformation
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/warper
    ${CMAKE_SOURCE_DIR}/toolboxes/image_io
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_image
  )

set(cpu_operator_header_files 
    cpuOperatorExport.h 
    hoIdentityOperator.h
    hoImageOperator.h
    hoDiagonalOperator.h
    hoDiagonalSumOperator.h
    hoFFTOperator.h
    hoPartialDerivativeOperator.h
    hoTvOperator.h
    hoTvPicsOperator.h 
    hoSPIRITOperator.h 
    hoSPIRIT2DOperator.h 
    hoSPIRIT2DTOperator.h 
    hoSPIRIT2DTDataFidelityOperator.h
    hoSPIRIT3DOperator.h 
    hoWaveletOperator.h 
    hoWavelet1DOperator.h 
    hoWavelet2DTOperator.h 
    hoWavelet3DOperator.h 
    hoMotionCompensation2DTOperator.h 
    ../curveFittingOperator.h 
    ../curveFittingCostFunction.h 
    twoParaExpDecayOperator.h 
    twoParaExpRecoveryOperator.h 
    threeParaExpRecoveryOperator.h 
    )

set(cpu_operator_src_files 
    hoSPIRITOperator.cpp 
    hoSPIRIT2DOperator.cpp 
    hoSPIRIT2DTOperator.cpp 
    hoSPIRIT2DTDataFidelityOperator.cpp 
    hoSPIRIT3DOperator.cpp 
    hoWaveletOperator.cpp 
    hoWavelet1DOperator.cpp 
    hoWavelet2DTOperator.cpp 
    hoWavelet3DOperator.cpp 
    hoMotionCompensation2DTOperator.cpp )

add_library(gadgetron_toolbox_cpuoperator SHARED ${operator_header_files} ${cpu_operator_header_files} ${cpu_operator_src_files} )

set_target_properties(gadgetron_toolbox_cpuoperator PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_cpuoperator 
                    gadgetron_toolbox_cpucore
                    gadgetron_toolbox_log
                    gadgetron_toolbox_cpucore_math 
                    gadgetron_toolbox_cpufft 
                    gadgetron_toolbox_cpuklt 
                    gadgetron_toolbox_cpudwt 
                    gadgetron_toolbox_mri_core 
                    gadgetron_toolbox_image_analyze_io )

install(TARGETS gadgetron_toolbox_cpuoperator DESTINATION lib COMPONENT main)

install(FILES 
  ${cpu_operator_header_files}
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
