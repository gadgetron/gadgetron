if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_GPUOPERATORS__)
endif ()



cuda_add_library(gadgetron_toolbox_gpuoperators SHARED 
  ../generalOperator.h
  ../linearOperator.h
  cuPartialDerivativeOperator.h
  cuLaplaceOperator.h
  cuTvOperator.h
  cuTv1dOperator.h
  cuConvolutionOperator.h
  cuPartialDerivativeOperator.cu
  cuPartialDerivativeOperator2.cu
  cuLaplaceOperator.cu
  cuTvOperator.cu
  cuTv1dOperator.cu
  cuConvolutionOperator.cu
  )

set_target_properties(gadgetron_toolbox_gpuoperators PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_gpuoperators
        gadgetron_toolbox_operator
  gadgetron_toolbox_gpucore 
  gadgetron_toolbox_gpunfft
  ${CUDA_LIBRARIES}
  ${CUDA_CUBLAS_LIBRARIES} 
  )


target_include_directories(gadgetron_toolbox_gpuoperators
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS gadgetron_toolbox_gpuoperators
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES 
  cuImageOperator.h
  cuDiagonalOperator.h
  cuDiagonalSumOperator.h
  cuPartialDerivativeOperator.h
  cuPartialDerivativeOperator2.h
  cuConvolutionOperator.h
  cuLaplaceOperator.h
  cuTvOperator.h
  cuTvPicsOperator.h
  cuTv1dOperator.h
  cuDownsampleOperator.h
  cuFFTOperator.h
  cuUpsampleOperator.h
  hoCuIdentityOperator.h
  hoCuPartialDerivativeOperator.h
  hoCuTvOperator.h
  hoCuTvPicsOperator.h
  hoCuEncodingOperatorContainer.h
  gpuoperators_export.h
  hoCuOperator.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)