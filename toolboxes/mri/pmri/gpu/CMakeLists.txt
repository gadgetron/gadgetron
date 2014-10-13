if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUPMRI__)
  ADD_DEFINITIONS(-DWIN32_LEAN_AND_MEAN)
endif (WIN32)

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories(  
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${CUDA_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
  ${ARMADILLO_INCLUDE_DIRS} 
  )

cuda_add_library(gadgetron_toolbox_gpuparallelmri SHARED 
    b1_map.h
    spirit_calibration.h
    cuCartesianSenseOperator.h
    cuNonCartesianKtSenseOperator.h
    cuNonCartesianSenseOperator.h
    cuSpiritOperator.h
    cuBuffer.h
    cuSenseBuffer.h
    cuSenseBufferCg.h
    cuSenseOperator.h
    gpupmri_export.h
    htgrappa.h
    senseOperator.h
    sense_utilities.h
    b1_map.cu
    b1_map_NIH_Souheil.cu
    spirit_calibration.cu
    sense_utilities.cu
    cuSenseOperator.cu
    cuCartesianSenseOperator.cu
    cuNonCartesianSenseOperator.cu
    cuNonCartesianKtSenseOperator.cu
    cuBuffer.cpp
    cuSenseBuffer.cpp
    cuSenseBufferCg.cpp
    cuSpiritBuffer.cpp
    htgrappa.cpp
    htgrappa.cu
  )

set_target_properties(gadgetron_toolbox_gpuparallelmri PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_gpuparallelmri 
  gadgetron_toolbox_gpucore 
  gadgetron_toolbox_gpunfft 
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_cpucore_math
  ${Boost_LIBRARIES}
  ${FFTW3_LIBRARIES} 
  ${CUDA_LIBRARIES} 
  ${CUDA_CUFFT_LIBRARIES} 
  ${CUDA_CUBLAS_LIBRARIES} 
  )

install(TARGETS gadgetron_toolbox_gpuparallelmri DESTINATION lib COMPONENT main)

install(FILES 
	b1_map.h
	sense_utilities.h
	htgrappa.h
	senseOperator.h
	cuSenseOperator.h
	cuCartesianSenseOperator.h
	cuNonCartesianSenseOperator.h
	cuNonCartesianKtSenseOperator.h
        cuSpiritOperator.h
        cuBuffer.h
	cuSenseBuffer.h
	cuSenseBufferCg.h
	cuSpiritBuffer.h
	gpupmri_export.h
DESTINATION include COMPONENT main)
