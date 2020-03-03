if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_CPUOPERATORS__)
endif ()



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

add_library(gadgetron_toolbox_cpuoperator SHARED ${operator_header_files} ${cpu_operator_header_files} ${cpu_operator_src_files} )

set_target_properties(gadgetron_toolbox_cpuoperator PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_include_directories(gadgetron_toolbox_cpuoperator
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(gadgetron_toolbox_cpuoperator
                    gadgetron_toolbox_operator
                    gadgetron_toolbox_cpucore
                    gadgetron_toolbox_cpureg
                    gadgetron_toolbox_log
                    gadgetron_toolbox_cpucore_math 
                    gadgetron_toolbox_cpufft 
                    gadgetron_toolbox_cpuklt 
                    gadgetron_toolbox_cpudwt
                    gadgetron_toolbox_mri_core
                    gadgetron_toolbox_image_analyze_io )

install(TARGETS gadgetron_toolbox_cpuoperator
		EXPORT gadgetron-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES 
  ${cpu_operator_header_files}
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
