
# matlab warpper
if (MATLAB_FOUND)

    message("MATLAB FOUND: matlab wrapper for gtplus toolbox will be compiled.")

    SET(CMAKE_DEBUG_POSTFIX)

    include_directories( 
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
                    optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
                    ${ISMRMRD_LIBRARIES} 
                    gadgetron_toolbox_matlab 
                    gadgetron_toolbox_image_analyze_io 
                    gadgetron_toolbox_cpucore 
                    gadgetron_toolbox_cpucore_math
                    gadgetron_toolbox_mri_core
                    gadgetron_toolbox_cpureg )

    if (WIN32)
        if ( HAS_64_BIT )
            SET(MATLAB_SUFFIX ".mexw64")
        else ( HAS_64_BIT )
            SET(MATLAB_SUFFIX ".mexw32")
        endif ( HAS_64_BIT )

        add_definitions(-DMATLAB_DLL_EXPORT_SYM)
    endif (WIN32)

    if (UNIX)
        if ( HAS_64_BIT )
            SET(MATLAB_SUFFIX ".mexa64")
        else ( HAS_64_BIT )
            SET(MATLAB_SUFFIX ".mexglx")
        endif ( HAS_64_BIT )
    endif (UNIX)

    if (APPLE)
        if ( HAS_64_BIT )
            SET(MATLAB_SUFFIX ".mexmaci64")
        endif ( HAS_64_BIT )
    endif(APPLE)

    if ( MKL_FOUND AND FFTW3_FOUND )

        add_library(Matlab_gt_read_analyze SHARED Matlab_gt_read_analyze.cpp)
        SET_TARGET_PROPERTIES(Matlab_gt_read_analyze PROPERTIES SUFFIX ${MATLAB_SUFFIX})
        install(TARGETS Matlab_gt_read_analyze DESTINATION ${GADGETRON_INSTALL_MATLAB_PATH} COMPONENT main)

        add_library(Matlab_gt_write_analyze SHARED Matlab_gt_write_analyze.cpp)
        SET_TARGET_PROPERTIES(Matlab_gt_write_analyze PROPERTIES SUFFIX ${MATLAB_SUFFIX})
        install(TARGETS Matlab_gt_write_analyze DESTINATION ${GADGETRON_INSTALL_MATLAB_PATH} COMPONENT main)

    endif ( MKL_FOUND AND FFTW3_FOUND )

else(MATLAB_FOUND)
    message("MATLAB NOT FOUND: matlab wrapper for gtplus toolbox will not be compiled.")
endif(MATLAB_FOUND)
