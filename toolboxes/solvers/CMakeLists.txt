

install(FILES 
  solver.h
  linearOperatorSolver.h
  cgSolver.h
  nlcgSolver.h
  lbfgsSolver.h
  lsqrSolver.h
  sbSolver.h
  sbcSolver.h
  cgCallback.h
  cgPreconditioner.h
  lwSolver.h
  lbfgsSolver.h
  gpSolver.h
  gpBbSolver.h
  eigenTester.h
  osMOMSolver.h
  osSPSSolver.h
  osLALMSolver.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

add_subdirectory(cpu)

if( CUDA_FOUND)
  add_subdirectory(gpu)
endif ()
