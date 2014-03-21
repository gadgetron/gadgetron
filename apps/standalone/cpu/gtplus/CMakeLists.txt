
# matlab warpper
if (MATLAB_FOUND)

    message("MATLAB FOUND: matlab wrapper for gtplus toolbox will be compiled.")

    SET(CMAKE_DEBUG_POSTFIX)

    include_directories( 
	  ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
	  ${MATLAB_INCLUDE_DIR}  
	  ${ISMRMRD_INCLUDE_DIR} 
	)

    link_directories(${Boost_LIBRARY_DIRS})
    link_libraries(${MATLAB_LIBRARIES} 
                    optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
                    gtplus 
                    cpucore 
                    cpucore_math)

    if (WIN32)
        if ( HAS_64_BIT )
            SET(MATLAB_SUFFIX ".mexw64")
        else ( HAS_64_BIT )
            SET(MATLAB_SUFFIX ".mexw32")
        endif ( HAS_64_BIT )
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

        # coil map 2D
        add_library(Matlab_compute_coil_map_2D SHARED Matlab_compute_coil_map_2D.cpp)
        SET_TARGET_PROPERTIES(Matlab_compute_coil_map_2D PROPERTIES SUFFIX ${MATLAB_SUFFIX})
        install(TARGETS Matlab_compute_coil_map_2D DESTINATION bin )

        # coil map 3D
        add_library(Matlab_compute_coil_map_3D SHARED Matlab_compute_coil_map_3D.cpp)
        SET_TARGET_PROPERTIES(Matlab_compute_coil_map_3D PROPERTIES SUFFIX ${MATLAB_SUFFIX})
        install(TARGETS Matlab_compute_coil_map_3D DESTINATION bin )

    endif ( MKL_FOUND AND FFTW3_FOUND)

else(MATLAB_FOUND)
    message("MATLAB NOT FOUND: matlab wrapper for gtplus toolbox will not be compiled.")
endif(MATLAB_FOUND)
