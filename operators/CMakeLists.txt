if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_SOLVERS__)
endif (WIN32)

install(FILES 	
	hoCuPartialDerivativeOperator.h
	subsetOperator.h
DESTINATION include)
