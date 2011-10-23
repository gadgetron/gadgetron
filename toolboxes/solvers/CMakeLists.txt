include_directories(		${ACE_INCLUDE_DIR}
					${CUDA_INCLUDE_DIRS}
					${CUDA_SDK_ROOT_DIR}/common/inc 
				       ${CMAKE_SOURCE_DIR}/toolboxes/ndarray
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpucore)

cuda_add_library(solvers SHARED 	cuCGPrecondWeights.cu
							cuPartialDerivativeOperator.cu
							cuLaplaceOperator.cu
							cuVariableGaussOperator.cu)

target_link_libraries(solvers gpucore ${CUDA_LIBRARIES}
${CUDA_CUBLAS_LIBRARIES} )

install(TARGETS solvers DESTINATION lib)

install(FILES 	
	solver.h
	cgSolver.h
	cuCGSolver.h
	sbSolver.h
	cuSBSolver.h
	cgPreconditioner.h
	cuCGPreconditioner.h
	cuCGPrecondWeights.h
	matrixOperator.h
	cuMatrixOperator_macros.h
	identityOperator.h
	cuIdentityOperator.h
	imageOperator.h
	cuImageOperator.h
	encodedImageOperator.h
	cuEncodedImageOperator.h
	partialDerivativeOperator.h
	cuPartialDerivativeOperator.h
	convolutionOperator.h
	cuConvolutionOperator.h
	laplaceOperator.h
	cuLaplaceOperator.h
	cuVariableGaussOperator.h
DESTINATION include)
