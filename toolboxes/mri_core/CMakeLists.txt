set(mri_core_header_files
        mri_core_def.h
        mri_core_utility.h
        mri_core_kspace_filter.h
        mri_core_grappa.h
        mri_core_spirit.h
        mri_core_coil_map_estimation.h
        mri_core_girf_correction.h
        mri_core_partial_fourier.h
        mri_core_stream.h
        )

set(mri_core_source_files
        mri_core_utility.cpp
        mri_core_grappa.cpp
        mri_core_spirit.cpp
        mri_core_kspace_filter.cpp
        mri_core_coil_map_estimation.cpp
        mri_core_girf_correction.cpp
        mri_core_partial_fourier.cpp
        mri_core_stream.cpp
        )

add_library(pingvin_toolbox_mri_core SHARED
        ${mri_core_header_files} ${mri_core_source_files})

set_target_properties(pingvin_toolbox_mri_core PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})
set_target_properties(pingvin_toolbox_mri_core PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(pingvin_toolbox_mri_core
        pingvin_core
        pingvin_toolbox_cpucore
        pingvin_toolbox_cpucore_math
        armadillo
        Boost::filesystem
        pingvin_toolbox_cpufft
        pingvin_toolbox_cpuklt
        pingvin_toolbox_image_analyze_io
        pingvin_toolbox_cpuklt)

target_include_directories(
        pingvin_toolbox_mri_core
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)


install(TARGETS pingvin_toolbox_mri_core
        EXPORT pingvin-export
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        COMPONENT main
        )


install(FILES
        ${mri_core_header_files}
        DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)

# if there is python, add python binding
if (BUILD_PYTHON_SUPPORT)
    add_library(pingvin_toolbox_mri_core_python MODULE
            mri_core_grappa_python.h
            mri_core_grappa_python.cpp
            mri_core_python.cpp)

    target_link_libraries(pingvin_toolbox_mri_core_python
            pingvin_toolbox_cpucore
            pingvin_toolbox_cpucore_math
            pingvin_toolbox_cpufft
            pingvin_toolbox_cpuklt
            pingvin_toolbox_image_analyze_io
            pingvin_toolbox_python
            pingvin_toolbox_mri_core
            )

    set_target_properties(pingvin_toolbox_mri_core_python PROPERTIES PREFIX "")

    install(TARGETS pingvin_toolbox_mri_core_python
            EXPORT pingvin-export
            DESTINATION ${PINGVIN_INSTALL_PYTHON_MODULE_PATH} COMPONENT main)
endif ()