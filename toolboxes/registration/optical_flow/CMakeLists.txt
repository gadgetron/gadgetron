install(FILES
  registrationSolver.h
  multiresRegistrationSolver.h
  opticalFlowSolver.h 
  resampleOperator.h
  opticalFlowOperator.h 
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

if(${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.819" )
  message("Compiling cpu based optical flow registration toolbox.")
  set(BUILD_CPU_OPTIMAL_FLOW_REG On)
else ()
  message("Armadillo (at least version 3.820) not found. Not compiling cpu-based optical flow registration toolbox. ")
  set(BUILD_CPU_OPTIMAL_FLOW_REG Off)
endif()

add_subdirectory(cpu)

if (CUDA_FOUND)
  message("Compiling gpu based optical flow registration toolbox.")
  add_subdirectory(gpu)
endif ()
