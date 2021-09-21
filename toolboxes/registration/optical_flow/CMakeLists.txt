install(FILES
  registrationSolver.h
  multiresRegistrationSolver.h
  opticalFlowSolver.h 
  resampleOperator.h
  opticalFlowOperator.h 
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

add_library(gadgetron_toolbox_registration_optical_flow INTERFACE)
target_link_libraries(gadgetron_toolbox_registration_optical_flow INTERFACE
        gadgetron_toolbox_cpu_solver)
target_include_directories(gadgetron_toolbox_registration_optical_flow
        INTERFACE
        $<INSTALL_INTERFACE:include/gadgetron>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

message("Compiling cpu based optical flow registration toolbox.")
set(BUILD_CPU_OPTIMAL_FLOW_REG On)

add_subdirectory(cpu)

if (CUDA_FOUND)
  message("Compiling gpu based optical flow registration toolbox.")
  add_subdirectory(gpu)
endif ()

install(TARGETS gadgetron_toolbox_registration_optical_flow
        EXPORT gadgetron-export
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        COMPONENT main
        )

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)