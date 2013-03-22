if (WIN32)
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
endif (WIN32)

include_directories( 
  ${Boost_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/toolboxes/cpucore 
  ${CMAKE_SOURCE_DIR}/toolboxes/cpucore/arma_math 
  ${CMAKE_SOURCE_DIR}/toolboxes/hostutils
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
  )

add_executable(cpu_denoise_TV denoise_TV.cpp)

target_link_libraries(cpu_denoise_TV cpucore cpucore_math hostutils)

install(TARGETS cpu_denoise_TV DESTINATION bin)
