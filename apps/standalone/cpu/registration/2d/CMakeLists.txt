include_directories( 
        ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/klt/cpu
        ${ACE_INCLUDE_DIR}
        ${ISMRMRD_INCLUDE_DIR}
    )

add_executable(register_HS_2d_cpu register_HS_2d.cpp)
add_executable(register_CK_2d_cpu register_CK_2d.cpp)

target_link_libraries(register_HS_2d_cpu 
  gadgetron_toolbox_hostutils 
  gadgetron_toolbox_cpureg 
  gadgetron_toolbox_cpucore 
  gadgetron_toolbox_cpucore_math
  armadillo
  )

target_link_libraries(register_CK_2d_cpu 
  gadgetron_toolbox_hostutils 
  gadgetron_toolbox_cpureg 
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_cpucore_math
  armadillo
  )

install(TARGETS 
  register_HS_2d_cpu
  register_CK_2d_cpu 
  DESTINATION bin COMPONENT main)
