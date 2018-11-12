if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_GPUNFFT__)
  add_definitions(-D_USE_MATH_DEFINES)
endif ()

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif()

include_directories( 
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_BINARY_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CUDA_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
  )

cuda_add_library(gadgetron_toolbox_gpunfft SHARED 
    cuNFFT.h
    gpunfft_export.h
    cuNFFT.cu
        cuNFFTOperator.cpp
  )

set_target_properties(gadgetron_toolbox_gpunfft PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_gpunfft 
  gadgetron_toolbox_gpufft
  gadgetron_toolbox_gpucore
  gadgetron_toolbox_log
  ${CUDA_LIBRARIES} 
  ${CUDA_CUFFT_LIBRARIES} 
  ${CUDA_CUBLAS_LIBRARIES}
  ${CUDA_CUSPARSE_LIBRARIES}
  )

install(TARGETS gadgetron_toolbox_gpunfft DESTINATION lib COMPONENT main)

install(FILES 
  cuNFFT.h
        ../NFFTOperator.h
  gpunfft_export.h 
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
