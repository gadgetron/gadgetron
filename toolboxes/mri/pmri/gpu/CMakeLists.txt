if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUPMRI__)
ADD_DEFINITIONS(-DWIN32_LEAN_AND_MEAN)
endif (WIN32)

if(WIN32)
link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

find_package(CULA REQUIRED)

include_directories(
  ${CUDA_INCLUDE_DIRS}
  ${CULA_INCLUDE_DIR} 
  ${Boost_INCLUDE_DIR} 
  #${CMAKE_SOURCE_DIR}/apps/gadgetron 
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  #${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/gadgets/core
  )

cuda_add_library(gpuparallelmri SHARED 	
  b1_map.cu
  sense_utilities.cu
  cuSenseOperator.cu
  cuCartesianSenseOperator.cu
  cuNonCartesianSenseOperator.cu
  cuNonCartesianKtSenseOperator.cu
  cuSenseRHSBuffer.cu
  htgrappa.cu
  )

target_link_libraries(gpuparallelmri 
  gpucore 
  gpunfft 
  ${Boost_LIBRARIES}
  ${FFTW3_LIBRARIES} 
  ${CUDA_LIBRARIES} 
  ${CUDA_CUFFT_LIBRARIES} 
  ${CUDA_CUBLAS_LIBRARIES} 
  ${CULA_LIBRARIES}
  )

install(TARGETS gpuparallelmri DESTINATION lib)

install(FILES 
	b1_map.h
	sense_utilities.h
	htgrappa.h
	senseOperator.h
	cuSenseOperator.h
	cuCartesianSenseOperator.h
	cuNonCartesianSenseOperator.h
	cuNonCartesianKtSenseOperator.h
	cuSenseRHSBuffer.h
	gpupmri_export.h
DESTINATION include)
