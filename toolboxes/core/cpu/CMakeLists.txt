if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUCORE__)
endif (WIN32)

include_directories(
  ${FFTW3_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  )

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

if (MKL_FOUND)
    INCLUDE_DIRECTORIES( ${MKL_INCLUDE_DIR} )
    LINK_DIRECTORIES( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
endif (MKL_FOUND)

add_library(cpucore ${LIBTYPE} 
  ../NDArray.h
  ../complext.h
  ../GadgetronException.h
  ../GadgetronCommon.h
  ../GadgetronTimer.h
  ../SerializableObject.h
  cpucore_export.h 
  hoNDArray.h
  hoNDArray.hxx
  hoNDArray_utils.h
  hoNDArray_fileio.h
  hoNDFFT.h
  ho2DArray.h
  ho2DArray.hxx
  ho3DArray.h
  ho3DArray.hxx
  ho4DArray.h
  ho4DArray.hxx
  ho5DArray.h
  ho5DArray.hxx
  ho6DArray.h
  ho6DArray.hxx
  ho7DArray.h
  ho7DArray.hxx 
  hoMatrix.h
  hoNDFFT.cpp
  )

target_link_libraries(cpucore 
  ${FFTW3_LIBRARIES} 
  ${Boost_LIBRARIES} 
  ${MKL_LIBRARIES} 
  )

install(TARGETS cpucore DESTINATION lib)

install(FILES
  cpucore_export.h 
  hoNDArray.h
  hoNDArray.hxx
  hoNDArray_utils.h
  hoNDArray_fileio.h
  hoNDFFT.h
  ho2DArray.h
  ho2DArray.hxx
  ho3DArray.h
  ho3DArray.hxx
  ho4DArray.h
  ho4DArray.hxx
  ho5DArray.h
  ho5DArray.hxx
  ho6DArray.h
  ho6DArray.hxx
  ho7DArray.h
  ho7DArray.hxx
  hoMatrix.h
  hoMatrix.cpp
  DESTINATION include)

if (ARMADILLO_FOUND)
    add_subdirectory(arma_math)
endif (ARMADILLO_FOUND)

add_subdirectory(hostutils)
