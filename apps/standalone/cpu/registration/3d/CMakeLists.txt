add_executable(register_CK_3d_cpu register_CK_3d.cpp)

target_link_libraries(register_CK_3d_cpu
  hostutils 
  cpureg 
  cpucore 
  cpucore_math
  ${ARMADILLO_LIBRARIES}
  )

install(TARGETS register_CK_3d_cpu DESTINATION bin)
