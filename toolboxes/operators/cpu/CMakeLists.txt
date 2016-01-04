if (WIN32)
    ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUOPERATORS__)
endif (WIN32)

include_directories(
    ${ISMRMRD_INCLUDE_DIR}
    ${FFTW3_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/operators
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/util
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/workflow
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/algorithm 
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/solver
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
    hoSPIRIT3DOperator.h )

set(cpu_operator_src_files 
    hoSPIRITOperator.cpp 
    hoSPIRIT2DOperator.cpp 
    hoSPIRIT2DTOperator.cpp 
    hoSPIRIT3DOperator.cpp )

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
                    gadgetron_toolbox_gtplus_io )

install(TARGETS gadgetron_toolbox_cpuoperator DESTINATION lib COMPONENT main)

install(FILES 
  ${cpu_operator_header_files}
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
