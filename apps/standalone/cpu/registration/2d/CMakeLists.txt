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
  
# matlab warpper
if (MATLAB_FOUND)
    message("MATLAB FOUND: ${MATLAB_INCLUDE_DIR}")
    if (WIN32)
    include_directories( ${MATLAB_INCLUDE_DIR} )
    add_library(Matlab_register_CK_2d_cpu SHARED Matlab_register_CK_2d.cpp)
    target_link_libraries(Matlab_register_CK_2d_cpu ${MATLAB_LIBRARIES} hostutils cpureg cpucore cpucore_math)
    if ( HAS_64_BIT )
        add_custom_command(TARGET Matlab_register_CK_2d_cpu 
                    POST_BUILD COMMAND del Matlab_register_CK_2d_cpu.mexw64 
                    POST_BUILD COMMAND rename Matlab_register_CK_2d_cpu.dll Matlab_register_CK_2d_cpu.mexw64 
                    WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}\\${CMAKE_CFG_INTDIR} )
    endif ( HAS_64_BIT )

    install(FILES Matlab_register_CK_2d_cpu.mexw64 DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )
  endif (WIN32)
else(MATLAB_FOUND)
  message("MATLAB NOT FOUND: no matlab warpper will be compiled ... ")
endif(MATLAB_FOUND)
