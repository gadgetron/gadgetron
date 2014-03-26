if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUFFT__)
endif (WIN32)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${FFTW3_INCLUDE_DIR}
  ${Boost_INCLUDE_DIR}
  )

add_library(cpufft ${LIBTYPE} 
  cpufft_export.h 
  hoNDFFT.h
  hoNDFFT.cpp
  )

target_link_libraries(cpufft
  cpucore
  ${FFTW3_LIBRARIES} 
  ${Boost_LIBRARIES} 
  )

install(TARGETS cpufft DESTINATION lib)

install(FILES
  cpufft_export.h 
  hoNDFFT.h
  DESTINATION include)
