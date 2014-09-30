if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUFFT__)
endif (WIN32)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${FFTW3_INCLUDE_DIR}
  ${Boost_INCLUDE_DIR}
  # ${ACE_INCLUDE_DIR}
  ${ISMRMRD_INCLUDE_DIR}
  ${ARMADILLO_INCLUDE_DIRS}
  )

add_library(cpufft ${LIBTYPE} 
  cpufft_export.h 
  hoNDFFT.h
  hoNDFFT.cpp
  )

target_link_libraries(cpufft
  cpucore
  cpucore_math
  ${FFTW3_LIBRARIES} 
  ${Boost_LIBRARIES} 
  )

install(TARGETS cpufft DESTINATION lib COMPONENT main)

install(FILES
  cpufft_export.h 
  hoNDFFT.h
  DESTINATION include COMPONENT main)
