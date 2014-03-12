if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_SOLVERS__)
endif (WIN32)



cuda_add_library( tomography_solvers SHARED
  cuNCGSolver.cu
)
install(FILES
	cuNCGSolver.h 	
	mlSolver.h
	cuMLSolver.h
	cuSARTSolver.h
	sartSolver.h
	nonLinearCGSolver.h
	hoGPBBSolver.h
DESTINATION include)
