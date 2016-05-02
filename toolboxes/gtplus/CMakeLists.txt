
if ( HAS_64_BIT )

    if (WIN32)
        ADD_DEFINITIONS(-D__BUILD_GADGETRON_PLUS__)
    endif (WIN32)

    if(WIN32)
        link_directories(${Boost_LIBRARY_DIRS})
    endif(WIN32)

    if (MKL_FOUND)
        list(APPEND EXTRA_MKL_LIBRARIES mkl_core)
        if ( USE_OPENMP )
            list(APPEND EXTRA_MKL_LIBRARIES mkl_intel_thread)
        endif ( USE_OPENMP )

        INCLUDE_DIRECTORIES( ${MKL_INCLUDE_DIR} )
        LINK_DIRECTORIES( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
    endif (MKL_FOUND)

    include_directories(
        ${ACE_INCLUDE_DIR} 
        ${Boost_INCLUDE_DIR}
        ${ISMRMRD_INCLUDE_DIR}
        ${FFTW3_INCLUDE_DIR}
        ${ARMADILLO_INCLUDE_DIRS}
        ${MKL_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/klt/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/dwt/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core
        ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
        ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
        ${CMAKE_SOURCE_DIR}/toolboxes/operators
        ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/solvers
        ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
        ${HDF5_INCLUDE_DIR}
        ${HDF5_INCLUDE_DIR}/cpp
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
        ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/gtplus
        ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
        ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
        ${CMAKE_SOURCE_DIR}/apps/gadgetron
        ${CMAKE_SOURCE_DIR}/apps/matlab
        ${CMAKE_SOURCE_DIR}/gadgets/mri_core 
    )

    set( gtplus_io_header_files 
        GtPlusIOExport.h 
        util/gtPlusIOBase.h
        util/gtPlusIOAnalyze.h)

    set( gtplus_io_src_files 
        util/gtPlusIOBase.cpp
        util/gtPlusIOAnalyze.cpp)

    add_library(gadgetron_toolbox_gtplus_io SHARED ${gtplus_io_header_files} ${gtplus_io_src_files} )
    set_target_properties(gadgetron_toolbox_gtplus_io PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
    target_link_libraries(gadgetron_toolbox_gtplus_io gadgetron_toolbox_log)

    install(TARGETS gadgetron_toolbox_gtplus_io DESTINATION lib COMPONENT main)

    install (FILES  
            ${gtplus_io_header_files} 
            DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

endif ( HAS_64_BIT )
