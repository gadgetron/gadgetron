if (WIN32)
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
endif (WIN32)

include_directories( 
  ${ACE_INCLUDE_DIR} 
  ${CMAKE_SOURCE_DIR}/toolboxes/core 
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu 
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math 
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
  )

add_executable(cpu_denoise_TV denoise_TV.cpp)

target_link_libraries(cpu_denoise_TV 
  cpucore 
  cpucore_math 
  hostutils
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
  )

install(TARGETS cpu_denoise_TV DESTINATION bin)
