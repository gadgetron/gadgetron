if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_MRI_CORE__)
endif ()


set(mri_core_header_files
        mri_core_export.h
        mri_core_def.h
        mri_core_data.h
        mri_core_utility.h
        mri_core_kspace_filter.h
        mri_core_grappa.h
        mri_core_spirit.h
        mri_core_coil_map_estimation.h
        mri_core_dependencies.h
        mri_core_acquisition_bucket.h
        mri_core_girf_correction.h
        mri_core_partial_fourier.h)

set(mri_core_source_files
        mri_core_utility.cpp
        mri_core_grappa.cpp
        mri_core_spirit.cpp
        mri_core_kspace_filter.cpp
        mri_core_coil_map_estimation.cpp
        mri_core_dependencies.cpp
        mri_core_girf_correction.cpp
        mri_core_partial_fourier.cpp)

add_library(gadgetron_toolbox_mri_core SHARED
        ${mri_core_header_files} ${mri_core_source_files})

set_target_properties(gadgetron_toolbox_mri_core PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
set_target_properties(gadgetron_toolbox_mri_core PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(gadgetron_toolbox_mri_core
        gadgetron_core
        gadgetron_toolbox_cpucore
        gadgetron_toolbox_cpucore_math
        ISMRMRD::ISMRMRD
        armadillo
        boost
        gadgetron_toolbox_cpufft
        gadgetron_toolbox_cpuklt
        gadgetron_toolbox_image_analyze_io
        gadgetron_toolbox_cpuklt)
target_include_directories(
        gadgetron_toolbox_mri_core
        PUBLIC
        $<INSTALL_INTERFACE:include/gadgetron>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)


install(TARGETS gadgetron_toolbox_mri_core
        EXPORT gadgetron-export
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        COMPONENT main
        )


install(FILES
        ${mri_core_header_files}
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

# if there is python, add python binding
if (Boost_PYTHON3_FOUND AND PYTHONLIBS_FOUND AND NUMPY_FOUND)


    add_library(gadgetron_toolbox_mri_core_python MODULE
            mri_core_grappa_python.h
            mri_core_grappa_python.cpp
            mri_core_python.cpp)

    target_link_libraries(gadgetron_toolbox_mri_core_python
            gadgetron_toolbox_cpucore
            gadgetron_toolbox_cpucore_math
            gadgetron_toolbox_cpufft
            gadgetron_toolbox_cpuklt
            gadgetron_toolbox_image_analyze_io
            gadgetron_toolbox_python
            gadgetron_toolbox_mri_core
            )

    if (WIN32)
        set_target_properties(gadgetron_toolbox_mri_core_python PROPERTIES SUFFIX .pyd)
    endif ()

    set_target_properties(gadgetron_toolbox_mri_core_python PROPERTIES PREFIX "")

    install(TARGETS gadgetron_toolbox_mri_core_python
            EXPORT gadgetron-export
            DESTINATION ${GADGETRON_INSTALL_PYTHON_MODULE_PATH} COMPONENT main)

endif ()

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)