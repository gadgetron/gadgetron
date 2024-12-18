install(FILES
  registrationSolver.h
  multiresRegistrationSolver.h
  opticalFlowSolver.h
  resampleOperator.h
  opticalFlowOperator.h
  DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)

add_library(pingvin_toolbox_registration_optical_flow INTERFACE)
target_link_libraries(pingvin_toolbox_registration_optical_flow INTERFACE
        pingvin_toolbox_cpu_solver)
target_include_directories(pingvin_toolbox_registration_optical_flow
        INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR})

message("Compiling cpu based optical flow registration toolbox.")
set(BUILD_CPU_OPTIMAL_FLOW_REG On)

add_subdirectory(cpu)