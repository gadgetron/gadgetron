if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUNFFT__)
  ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
endif (WIN32)

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories( 
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CUDA_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
  )

cuda_add_library(gadgetron_toolbox_gpunfft SHARED 
    cuNFFT.h
    cuNFFTOperator.h
    gpunfft_export.h
    cuNFFT.cu 
    cuNFFTOperator.cu
  )

set_target_properties(gadgetron_toolbox_gpunfft PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_gpunfft 
  gadgetron_toolbox_gpufft
  gadgetron_toolbox_gpucore 
  ${Boost_LIBRARIES}
  ${CUDA_LIBRARIES} 
  ${CUDA_CUFFT_LIBRARIES} 
  ${CUDA_CUBLAS_LIBRARIES}
  )

install(TARGETS gadgetron_toolbox_gpunfft DESTINATION lib COMPONENT main)

install(FILES 
  cuNFFT.h 
  cuNFFTOperator.h 
  gpunfft_export.h 
  DESTINATION include COMPONENT main)
