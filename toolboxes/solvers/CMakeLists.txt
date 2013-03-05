if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_SOLVERS__)
endif (WIN32)

include_directories(		${CUDA_INCLUDE_DIRS}
				${Boost_INCLUDE_DIR}
				${CMAKE_SOURCE_DIR}/toolboxes/cpucore
				${CMAKE_SOURCE_DIR}/toolboxes/operators
				${CMAKE_SOURCE_DIR}/toolboxes/solvers
				${CMAKE_SOURCE_DIR}/toolboxes/gpucore)

cuda_add_library(solvers SHARED cuCgPrecondWeights.cu cuGPBBSolver.cu)

target_link_libraries(solvers gpucore ${CUDA_LIBRARIES}
${CUDA_CUBLAS_LIBRARIES} )

install(TARGETS solvers DESTINATION lib)

install(FILES 	
	solver.h
	linearSolver.h
	cgSolver.h
	cuCgSolver.h
	sbSolver.h
	sbcSolver.h
	cuSbSolver_macros.h
	cuSbCgSolver.h
	cuSbcCgSolver.h
	cgCallback.h	
	cgPreconditioner.h
	cuCgPreconditioner.h
	cuCgPrecondWeights.h
	lwSolver.h
	cuLwSolver.h
	cuSbLwSolver.h
	cuSbcLwSolver.h
	linearOperator.h
	cuGPBBSolver.h
	gpsolver.h
	gpBBSolver.h
	solvers_export.h
DESTINATION include)
