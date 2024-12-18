add_library(pingvin_toolbox_cpunfft SHARED
    hoNFFT.h
    hoNFFT.cpp
    ConvolutionMatrix.h
    ConvolutionMatrix.cpp
    hoGriddingConvolution.h
    hoGriddingConvolution.cpp
	  hoNFFTOperator.cpp
)

set_target_properties(pingvin_toolbox_cpunfft PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_include_directories(pingvin_toolbox_cpunfft
        PUBLIC
		$<INSTALL_INTERFACE:include>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)

target_link_libraries(pingvin_toolbox_cpunfft
    pingvin_toolbox_cpufft
    pingvin_toolbox_cpucore
    pingvin_toolbox_log
	pingvin_toolbox_operator
  )

install(TARGETS pingvin_toolbox_cpunfft
		EXPORT pingvin-export
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
		RUNTIME DESTINATION bin
		COMPONENT main
		)

install(FILES
    hoNFFT.h
    ConvolutionMatrix.h
    hoGriddingConvolution.h
    DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)