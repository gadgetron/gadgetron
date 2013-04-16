if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUCORE_MATH__)
endif (WIN32)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${ARMADILLO_INCLUDE_DIRS}
  )

add_library(cpucore_math SHARED 
  hoNDArray_operators.cpp
  hoNDArray_elemwise.cpp
  hoNDArray_blas.cpp
  )

target_link_libraries(cpucore_math 
  ${ARMADILLO_LIBRARIES}
  )

install(TARGETS cpucore_math DESTINATION lib)

install(FILES 	
  cpucore_math_export.h
  hoArmadillo.h
  hoNDArray_operators.h
  hoNDArray_elemwise.h
  hoNDArray_blas.h
  DESTINATION include)
