include_directories(   
  ${Boost_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
)

install(FILES
  registrationSolver.h
  multiresRegistrationSolver.h
  resampleOperator.h
  DESTINATION include)

add_subdirectory(cpu)

if (CUDA_FOUND)
  add_subdirectory(gpu)
endif (CUDA_FOUND)
