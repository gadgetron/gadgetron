if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_GPUFFT__)
endif ()

include_directories( 
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CUDA_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
)

cuda_add_library(gadgetron_toolbox_gpufft SHARED 
    cuNDFFT.h
    cuNDFFT.cpp
    cuNDFFT.cu
  )

set_target_properties(gadgetron_toolbox_gpufft PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_gpufft 
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_log
  gadgetron_toolbox_gpucore 
  ${CUDA_LIBRARIES} 
  ${CUDA_CUFFT_LIBRARIES} 
  )

install(TARGETS gadgetron_toolbox_gpufft DESTINATION lib COMPONENT main)

install(FILES
  gpufft_export.h
  cuNDFFT.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
