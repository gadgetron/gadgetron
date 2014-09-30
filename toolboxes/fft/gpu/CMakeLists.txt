if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUFFT__)
endif (WIN32)

include_directories( 
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CUDA_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
)

cuda_add_library(gpufft SHARED 
    cuNDFFT.h
    cuNDFFT.cpp
  )

target_link_libraries(gpufft 
  cpucore
  gpucore 
  ${Boost_LIBRARIES}
  ${CUDA_LIBRARIES} 
  ${CUDA_CUFFT_LIBRARIES} 
  )

install(TARGETS gpufft DESTINATION lib COMPONENT main)

install(FILES
  gpufft_export.h
  cuNDFFT.h
  DESTINATION include COMPONENT main)
