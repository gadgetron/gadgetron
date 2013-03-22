include_directories(
  ${Boost_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/cpucore
  )

install(FILES 	
  solver.h
  linearOperatorSolver.h
  cgSolver.h
  sbSolver.h
  sbcSolver.h
  cgCallback.h	
  cgPreconditioner.h
  lwSolver.h
  gpSolver.h
  gpBBSolver.h
  DESTINATION include)

add_subdirectory(cpu)
#add_subdirectory(gpu)
