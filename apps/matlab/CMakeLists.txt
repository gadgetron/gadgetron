# matlab warpper
if (MATLAB_FOUND)

    message("MATLAB FOUND: matlab wrapper for toolbox will be compiled.")

    set(CMAKE_DEBUG_POSTFIX)

    include_directories( 
        ${Boost_INCLUDE_DIR} 
        ${ARMADILLO_INCLUDE_DIRS} 
        ${ACE_INCLUDE_DIR} 
        ${ISMRMRD_INCLUDE_DIR}
        ${FFTW3_INCLUDE_DIR}
        ${MKL_INCLUDE_DIR}
        ${MATLAB_INCLUDE_DIR}
        ${CMAKE_BINARY_DIR}/toolboxes/core
        ${CMAKE_SOURCE_DIR}/toolboxes/matlab
        ${CMAKE_SOURCE_DIR}/toolboxes/core
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
        ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/klt/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/dwt/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/ffd
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
        ${CMAKE_SOURCE_DIR}/toolboxes/mri_image
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/application
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/dissimilarity
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/register
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/solver
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/transformation
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/warper
        ${CMAKE_SOURCE_DIR}/toolboxes/operators
        ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/solvers
        ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
        ${CMAKE_SOURCE_DIR}/toolboxes/image_io
        ${CMAKE_SOURCE_DIR}/toolboxes/image
        ${CMAKE_SOURCE_DIR}/toolboxes/image/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/cmr
    )

    link_directories(${Boost_LIBRARY_DIRS})
    link_libraries(${MATLAB_LIBRARIES} 
                    optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
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
                    gadgetron_toolbox_cmr 
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


    add_library(Matlab_gt_deformation_field_reg_3D SHARED matlab_gt_deformation_field_reg_3D.cpp)
    SET_TARGET_PROPERTIES(Matlab_gt_deformation_field_reg_3D PROPERTIES SUFFIX ${MATLAB_SUFFIX} PREFIX "")
    install(TARGETS Matlab_gt_deformation_field_reg_3D DESTINATION ${GADGETRON_INSTALL_MATLAB_PATH} )

else()
    message("MATLAB NOT FOUND: matlab wrapper for toolbox will not be compiled.")
endif()