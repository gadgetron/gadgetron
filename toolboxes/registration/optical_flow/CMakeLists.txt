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

if(ARMADILLO_FOUND)
  message("Compiling optical flow registration toolbox.")
  add_subdirectory(cpu)
else (ARMADILLO_FOUND)
  message("Armadillo not found. Not compiling the optical flow registration toolbox. ")  
endif(ARMADILLO_FOUND)

if (CUDA_FOUND)
  add_subdirectory(gpu)
endif (CUDA_FOUND)
