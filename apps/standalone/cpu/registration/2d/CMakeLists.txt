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

# matlab wrapper
if (MATLAB_FOUND)
  message("Matlab found> ${MATLAB_INCLUDE_DIR}. Matlab registration wrapper is being compiled.")
  set(CMAKE_DEBUG_POSTFIX)
  if (WIN32)
    include_directories( ${MATLAB_INCLUDE_DIR} )

    if (WIN32)
        if ( HAS_64_BIT )
            set(MATLAB_SUFFIX ".mexw64")
        else ()
            set(MATLAB_SUFFIX ".mexw32")
        endif ()

        add_definitions(-DMATLAB_DLL_EXPORT_SYM)
    endif ()

    if (UNIX)
        if ( HAS_64_BIT )
            set(MATLAB_SUFFIX ".mexa64")
        else ()
            set(MATLAB_SUFFIX ".mexglx")
        endif ()
    endif ()

    if (APPLE)
        if ( HAS_64_BIT )
            set(MATLAB_SUFFIX ".mexmaci64")
        endif ()
    endif()
    
    add_library(Matlab_register_CK_2d_cpu SHARED Matlab_register_CK_2d.cpp)
    target_link_libraries(Matlab_register_CK_2d_cpu ${MATLAB_LIBRARIES} gadgetron_toolbox_hostutils gadgetron_toolbox_cpureg gadgetron_toolbox_cpucore gadgetron_toolbox_cpucore_math)
    if ( HAS_64_BIT )
      set_target_properties(Matlab_register_CK_2d_cpu PROPERTIES SUFFIX ${MATLAB_SUFFIX})
    endif ()    
    install(TARGETS Matlab_register_CK_2d_cpu DESTINATION ${GADGETRON_INSTALL_MATLAB_PATH} COMPONENT main)
  endif ()
else()
  message("Matlab not found. Matlab wrapper for registration toolbox will not be compiled.")
endif()
