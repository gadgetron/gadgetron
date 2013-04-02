if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUSOLVERS__)
endif (WIN32)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/cpucore/
  ${CMAKE_SOURCE_DIR}/toolboxes/cpucore/arma_math
  )

install(FILES 	
  hoCgSolver.h
  hoSbCgSolver.h
  hoGPBBSolver.h
  hoCgPreconditioner.h
  DESTINATION include)
