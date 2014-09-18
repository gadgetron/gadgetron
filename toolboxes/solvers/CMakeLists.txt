include_directories(
  ${Boost_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  )

install(FILES 	
  solver.h
  linearOperatorSolver.h
  cgSolver.h
  nlcgSolver.h
  sbSolver.h
  sbcSolver.h
  cgCallback.h	
  cgPreconditioner.h
  lwSolver.h
  lbfgsSolver.h
  gpSolver.h
  gpBbSolver.h
  DESTINATION include)

IF(ARMADILLO_FOUND)
  add_subdirectory(cpu)
ENDIF(ARMADILLO_FOUND)

IF( CUDA_FOUND)
  add_subdirectory(gpu)
ENDIF (CUDA_FOUND)
