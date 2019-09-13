
# matlab warpper
if (MATLAB_FOUND)

    message("MATLAB FOUND: matlab wrapper for gtplus toolbox will be compiled.")

    set(CMAKE_DEBUG_POSTFIX)

    include_directories( 
        ${CMAKE_BINARY_DIR}/toolboxes/core
        ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/klt/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
        ${CMAKE_SOURCE_DIR}/toolboxes/matlab
        ${CMAKE_SOURCE_DIR}/toolboxes/image_io
        ${ACE_INCLUDE_DIR}
        ${MATLAB_INCLUDE_DIR}  
        ${ISMRMRD_INCLUDE_DIR} 
    )

    link_directories(${Boost_LIBRARY_DIRS})
    link_libraries(${MATLAB_LIBRARIES}
                    ${ISMRMRD_LIBRARIES}
                    ${FFTW3_LIBRARIES}
                    ${Boost_LIBRARIES}
                    ${ARMADILLO_LIBRARIES}
                    gadgetron_toolbox_matlab 
                    gadgetron_toolbox_cpucore 
                    gadgetron_toolbox_cpucore_math
                    gadgetron_toolbox_mri_core
                    gadgetron_toolbox_cpureg 
                    gadgetron_toolbox_cpufft 
                    gadgetron_toolbox_cpuklt 
                    gadgetron_toolbox_cpuoperator 
                    gadgetron_toolbox_log )

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

    add_library(Matlab_gt_read_analyze SHARED Matlab_gt_read_analyze.cpp)
    set_target_properties(Matlab_gt_read_analyze PROPERTIES SUFFIX ${MATLAB_SUFFIX} PREFIX "")
    install(TARGETS Matlab_gt_read_analyze DESTINATION ${GADGETRON_INSTALL_MATLAB_PATH} COMPONENT main)

    add_library(Matlab_gt_write_analyze SHARED Matlab_gt_write_analyze.cpp)
    set_target_properties(Matlab_gt_write_analyze PROPERTIES SUFFIX ${MATLAB_SUFFIX} PREFIX "")
    install(TARGETS Matlab_gt_write_analyze DESTINATION ${GADGETRON_INSTALL_MATLAB_PATH} COMPONENT main)

else()
    message("MATLAB NOT FOUND: matlab wrapper for gtplus toolbox will not be compiled.")
endif()
