if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_SOLVERS__)
endif (WIN32)

include_directories(		${CUDA_INCLUDE_DIRS}
				${Boost_INCLUDE_DIR}
				${CMAKE_SOURCE_DIR}/toolboxes/cpucore
				${CMAKE_SOURCE_DIR}/toolboxes/solvers
				${CMAKE_SOURCE_DIR}/toolboxes/gpucore)

cuda_add_library(operators SHARED cuPartialDerivativeOperator.cu
                                cuLaplaceOperator.cu
                                cuVariableGaussOperator.cu)

target_link_libraries(operators gpucore ${CUDA_LIBRARIES}
${CUDA_CUBLAS_LIBRARIES} )

install(TARGETS operators DESTINATION lib)

install(FILES 	
	linearOperator.h
	cuLinearOperator_macros.h
	identityOperator.h
	imageOperator.h
	cuImageOperator.h
	diagonalOperator.h
	cuDiagonalOperator.h
	encodedImageOperator.h
	cuEncodedImageOperator.h
	partialDerivativeOperator.h
	cuPartialDerivativeOperator.h
	convolutionOperator.h
	cuConvolutionOperator.h
	laplaceOperator.h
	cuLaplaceOperator.h
	cuVariableGaussOperator.h
	encodingOperatorContainer.h
DESTINATION include)
