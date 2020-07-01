if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_GPUNFFT__)
  add_definitions(-D_USE_MATH_DEFINES)
endif ()



cuda_add_library(gadgetron_toolbox_gpunfft SHARED 
    cuNFFT.h
    gpunfft_export.h
    cuNFFT.cu 
    cuNFFTOperator.cpp
  )

target_include_directories(gadgetron_toolbox_gpunfft
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(gadgetron_toolbox_gpunfft PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_gpunfft
  gadgetron_toolbox_nfft
  gadgetron_toolbox_operator
  gadgetron_toolbox_gpufft
  gadgetron_toolbox_gpucore
  gadgetron_toolbox_log
  ${CUDA_LIBRARIES} 
  ${CUDA_CUFFT_LIBRARIES} 
  ${CUDA_CUBLAS_LIBRARIES}
  ${CUDA_CUSPARSE_LIBRARIES}
  )
install(TARGETS gadgetron_toolbox_gpunfft
		EXPORT gadgetron-export
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
		RUNTIME DESTINATION bin
		COMPONENT main
		)


install(FILES 
  cuNFFT.h 
  gpunfft_export.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)