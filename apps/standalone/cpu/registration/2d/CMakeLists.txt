include_directories( 
        ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
        ${ACE_INCLUDE_DIR}
        ${ISMRMRD_INCLUDE_DIR}
    )

add_executable(register_HS_2d_cpu register_HS_2d.cpp)
add_executable(register_CK_2d_cpu register_CK_2d.cpp)

target_link_libraries(register_HS_2d_cpu 
  hostutils 
  cpureg 
  cpucore 
  cpucore_math
  ${ARMADILLO_LIBRARIES}
  )

target_link_libraries(register_CK_2d_cpu 
  hostutils 
  cpureg 
  cpucore
  cpucore_math
  ${ARMADILLO_LIBRARIES}
  )

install(TARGETS 
  register_HS_2d_cpu
  register_CK_2d_cpu 
  DESTINATION bin)

# matlab wrapper
if (MATLAB_FOUND)
  message("Matlab found> ${MATLAB_INCLUDE_DIR}. Matlab registration wrapper is being compiled.")
  SET(CMAKE_DEBUG_POSTFIX)
  if (WIN32)
    include_directories( ${MATLAB_INCLUDE_DIR} )
    add_library(Matlab_register_CK_2d_cpu SHARED Matlab_register_CK_2d.cpp)
    target_link_libraries(Matlab_register_CK_2d_cpu ${MATLAB_LIBRARIES} hostutils cpureg cpucore cpucore_math)
    if ( HAS_64_BIT )				
      SET_TARGET_PROPERTIES(Matlab_register_CK_2d_cpu PROPERTIES SUFFIX .mexw64)
    endif ( HAS_64_BIT )    
    install(TARGETS Matlab_register_CK_2d_cpu DESTINATION bin )
  endif (WIN32)
else(MATLAB_FOUND)
  message("Matlab not found. Matlab wrapper for registration toolbox will not be compiled.")
endif(MATLAB_FOUND)
