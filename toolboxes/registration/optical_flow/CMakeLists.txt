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
  opticalFlowSolver.h 
  resampleOperator.h
  opticalFlowOperator.h 
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

if(ARMADILLO_FOUND)
  if(${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.819" )
    message("Compiling cpu based optical flow registration toolbox.")
    set(BUILD_CPU_OPTIMAL_FLOW_REG On)
  else (${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.819" )
    message("Armadillo (at least version 3.820) not found. Not compiling cpu-based optical flow registration toolbox. ")  
    set(BUILD_CPU_OPTIMAL_FLOW_REG Off)
  endif(${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.819" )
endif(ARMADILLO_FOUND)

add_subdirectory(cpu)

if (CUDA_FOUND)
  message("Compiling gpu based optical flow registration toolbox.")
  add_subdirectory(gpu)
endif (CUDA_FOUND)
