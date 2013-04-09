include_directories(   
  ${Boost_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
)

install(FILES
  registrationSolver.h
  multiresRegistrationSolver.h
  resampleOperator.h
  DESTINATION include)

if(ARMADILLO_FOUND AND ${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.39" )
  add_subdirectory(cpu)
endif(ARMADILLO_FOUND AND ${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.39" )

if (CUDA_FOUND)
  add_subdirectory(gpu)
endif (CUDA_FOUND)
