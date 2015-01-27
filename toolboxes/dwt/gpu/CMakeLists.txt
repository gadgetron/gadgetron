if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUDWT__)
endif (WIN32)

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)


include_directories( 
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CUDA_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
)

cuda_add_library(gadgetron_toolbox_gpudwt SHARED 
    cuNDDWT.h
    cuNDDWT.cu
  )

set_target_properties(gadgetron_toolbox_gpudwt PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_gpudwt 
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_gpucore 
  ${Boost_LIBRARIES}
  ${CUDA_LIBRARIES} 
  ${CUDA_CUFFT_LIBRARIES} 
  )

install(TARGETS gadgetron_toolbox_gpudwt DESTINATION lib COMPONENT main)

install(FILES
  cuNDDWT.h
  DESTINATION include COMPONENT main)
