if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_SOLVERS__)
endif (WIN32)


cuda_add_library( tomography_operators SHARED
  cuWaveletOperator.cu
)

add_executable(waveletTest waveletTest.cpp )


target_link_libraries(waveletTest tomography_operators gpucore cpucore_math gpuoperators gpusolvers hostutils ${HDF5_LIBRARIES} ${CUDA_LIBRARIES} ${CULA_LIBRARIES} ${Boost_LIBRARIES})


install(FILES 	
	hoCuPartialDerivativeOperator.h
	subsetOperator.h
DESTINATION include)
