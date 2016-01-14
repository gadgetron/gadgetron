if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_TOOLBOX_HYPER__)
endif (WIN32)

if(WIN32)
link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories(
  ${CUDA_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CUDA_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
  )

cuda_add_library(gadgetron_toolbox_hyper SHARED 
    CSIOperator.cpp
    CSI_utils.cu
	CSIOperator.h
	CSI_utils.h
	gadgetron_toolbox_hyper_export.h
  )

set_target_properties(gadgetron_toolbox_hyper PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_hyper 
  gadgetron_toolbox_gpucore 
  gadgetron_toolbox_gpunfft
  ${Boost_LIBRARIES}
  ${CUDA_CUDART_LIBRARY}
  ${CUDA_CUBLAS_LIBRARIES} 
  )

install(TARGETS gadgetron_toolbox_hyper DESTINATION lib COMPONENT main)

install(FILES 
    CSIOperator.h
    CSI_utils.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
