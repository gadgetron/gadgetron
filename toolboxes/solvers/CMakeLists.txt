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
  gpSolver.h
  gpBbSolver.h
  eigenTester.h
  osMOMSolver.h
  osSPSSolver.h
  osLALMSolver.h

  DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)

add_library(pingvin_toolbox_solvers INTERFACE)
target_include_directories(pingvin_toolbox_solvers
        INTERFACE
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
install(TARGETS pingvin_toolbox_solvers EXPORT pingvin-export)

add_subdirectory(cpu)
