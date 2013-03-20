if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUCORE_MATH__)
endif (WIN32)

include_directories(
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
  hoNDArray_fileio.h
  hoNDArray_elemwise.h
  hoNDArray_blas.h
  hoArmadillo.h
  DESTINATION include)
