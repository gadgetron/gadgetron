add_executable(register_CK_3d_cpu register_CK_3d.cpp)

target_link_libraries(register_CK_3d_cpu
  gadgetron_toolbox_hostutils 
  gadgetron_toolbox_cpureg 
  gadgetron_toolbox_cpucore 
  gadgetron_toolbox_cpucore_math
  ${ARMADILLO_LIBRARIES}
  )

install(TARGETS register_CK_3d_cpu DESTINATION bin COMPONENT main)
