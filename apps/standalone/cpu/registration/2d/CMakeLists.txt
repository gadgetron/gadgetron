add_executable(register_HS_2d_cpu register_HS_2d.cpp)
add_executable(register_CK_2d_cpu register_CK_2d.cpp)

target_link_libraries(register_HS_2d_cpu 
  hostutils 
  cpureg 
  cpucore 
  cpucore_math
  #cpuoperators
  #cpusolvers 
  )

target_link_libraries(register_CK_2d_cpu 
  hostutils 
  cpureg 
  cpucore
  cpucore_math
#  cpuoperators
#  cpusolvers 
  )

install(TARGETS 
  register_HS_2d_cpu
  register_CK_2d_cpu 
  DESTINATION bin)
