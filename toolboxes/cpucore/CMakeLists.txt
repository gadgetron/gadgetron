if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUCORE__)
endif (WIN32)

include_directories(
  ${ACE_INCLUDE_DIR} 
  ${Boost_INCLUDE_DIR}
  ${FFTW3_INCLUDE_DIR}
  )

configure_file(cpucore_defines.h.in ${CMAKE_CURRENT_SOURCE_DIR}/cpucore_defines.h)

add_library(cpucore SHARED 
  FFT.cpp
  )

target_link_libraries(cpucore 
  ${FFTW3_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
  )

install(TARGETS cpucore DESTINATION lib)

install(FILES 	
  cpucore_export.h
  cpucore_defines.h
  NDArray.h
  hoNDArray.h
  hoNDArray_fileio.h
  complext.h
  vector_td.h
  vector_td_operators.h
  vector_td_utilities.h
  real_utilities.h
  GadgetronException.h
  GadgetronTimer.h
  FFT.h
  DESTINATION include)

if (ARMADILLO_FOUND)
add_subdirectory(math)
endif (ARMADILLO_FOUND)
