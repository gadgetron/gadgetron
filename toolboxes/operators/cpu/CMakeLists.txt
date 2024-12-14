set(cpu_operator_header_files
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
    ../curveFittingOperator.h
    ../curveFittingCostFunction.h
    twoParaExpDecayOperator.h
    twoParaExpRecoveryOperator.h
    threeParaExpRecoveryOperator.h
    hoMotionCompensation2DTOperator.h
    )

set(cpu_operator_src_files
    hoSPIRITOperator.cpp
    hoSPIRIT2DOperator.cpp
    hoSPIRIT2DTOperator.cpp
    hoSPIRIT2DTDataFidelityOperator.cpp
    hoMotionCompensation2DTOperator.cpp
    hoSPIRIT3DOperator.cpp
    hoWaveletOperator.cpp
    hoWavelet1DOperator.cpp
    hoWavelet2DTOperator.cpp
    hoWavelet3DOperator.cpp
    )

add_library(pingvin_toolbox_cpuoperator SHARED ${operator_header_files} ${cpu_operator_header_files} ${cpu_operator_src_files} )

set_target_properties(pingvin_toolbox_cpuoperator PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_include_directories(pingvin_toolbox_cpuoperator
        PUBLIC
        $<INSTALL_INTERFACE:${PINGVIN_INSTALL_INCLUDE_PATH}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(pingvin_toolbox_cpuoperator
                    pingvin_toolbox_operator
                    pingvin_toolbox_cpucore
                    pingvin_toolbox_cpureg
                    pingvin_toolbox_log
                    pingvin_toolbox_cpucore_math
                    pingvin_toolbox_cpufft
                    pingvin_toolbox_cpuklt
                    pingvin_toolbox_cpudwt
                    pingvin_toolbox_mri_core
                    pingvin_toolbox_image_analyze_io )

install(TARGETS pingvin_toolbox_cpuoperator
		EXPORT pingvin-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
  ${cpu_operator_header_files}
  DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)