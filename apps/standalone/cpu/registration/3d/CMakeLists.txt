add_executable(register_CK_3d_cpu register_CK_3d.cpp)

target_link_libraries(register_CK_3d_cpu
  hostutils 
  cpureg 
  cpucore 
  cpucore_math
#  cpuoperators
#  cpusolvers 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
  )

install(TARGETS register_CK_3d_cpu DESTINATION bin)
