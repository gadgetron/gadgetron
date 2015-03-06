if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUSOLVERS__)
endif (WIN32)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
  )

install(FILES 	
  hoCgSolver.h
  hoSbCgSolver.h
  hoGpBbSolver.h
  hoCgPreconditioner.h
  hoSolverUtils.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
