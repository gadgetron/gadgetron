if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUCORE_MATH__)
endif (WIN32)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${ARMADILLO_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}/cpp
  )

if (MKL_FOUND)
    INCLUDE_DIRECTORIES( ${MKL_INCLUDE_DIR} )
    LINK_DIRECTORIES( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
endif (MKL_FOUND)

add_library(cpucore_math ${LIBTYPE} 
    hoNDArray_operators.h
    hoNDArray_operators.cpp
    hoNDArray_elemwise.h
    hoNDArray_elemwise.cpp
    hoNDArray_blas.h
    hoNDArray_blas.cpp
    hoNDArray_reductions.cpp
    )

target_link_libraries(cpucore_math 
    cpucore
    ${ARMADILLO_LIBRARIES}
    ${MKL_LIBRARIES}
    )

install(TARGETS cpucore_math DESTINATION lib)

install(FILES 	
  cpucore_math_export.h
  hoArmadillo.h
  hoNDArray_operators.h
  hoNDArray_elemwise.h
  hoNDArray_blas.h
  hoNDArray_reductions.h
  hoNDArray_math.h
  DESTINATION include)
