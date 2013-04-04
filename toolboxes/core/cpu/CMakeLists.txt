if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUCORE__)
endif (WIN32)

include_directories(
  ${FFTW3_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  )

add_library(cpucore SHARED 
  hoFFT.cpp
  )

target_link_libraries(cpucore 
  ${FFTW3_LIBRARIES} ${Boost_LIBRARIES}
  )

install(TARGETS cpucore DESTINATION lib)

install(FILES 	
  hoNDArray.h
  hoNDArray_utils.h
  hoNDArray_fileio.h
  hoFFT.h
  DESTINATION include)

if (ARMADILLO_FOUND)
add_subdirectory(arma_math)
endif (ARMADILLO_FOUND)

add_subdirectory(hostutils)
