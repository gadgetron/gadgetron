

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

add_library(gadgetron_toolbox_solvers INTERFACE)
target_include_directories(gadgetron_toolbox_solvers
        INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(cpu)

if( CUDA_FOUND)
  add_subdirectory(gpu)
endif ()

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)