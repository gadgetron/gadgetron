if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_DEBUGGING__)
endif ()

include_directories(${HDF5_INCLUDE_DIRS})

set(gadgetron_debugging_header_files
        AcquisitionPassthroughGadget.h
        CplxDumpGadget.h
        ImageWriterGadget.h
        PseudoReplicatorGadget.h
        RateLimitGadget.h
        WhiteNoiseInjectorGadget.h
        )

set(gadgetron_debugging_src_files
        AcquisitionPassthroughGadget.cpp
        CplxDumpGadget.cpp
        ImageWriterGadget.cpp
        PseudoReplicatorGadget.cpp
        RateLimitGadget.cpp
        WhiteNoiseInjectorGadget.cpp
        )

set(gadgetron_debugging_config_files
        )

add_library(gadgetron_debugging SHARED
    gadgetron_debugging_export.h
    ${gadgetron_debugging_header_files}
    ${gadgetron_debugging_src_files}
    ${gadgetron_debuggingconfig_files}
    )

set_target_properties(gadgetron_debugging PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_debugging
        gadgetron_core
        gadgetron_mricore
        gadgetron_core_writers
        gadgetron_toolbox_log
        gadgetron_toolbox_cpucore
        gadgetron_toolbox_cpufft
        gadgetron_toolbox_image_analyze_io
        gadgetron_toolbox_denoise
        gadgetron_toolbox_hostutils
        gadgetron_toolbox_cpuoperator
        gadgetron_toolbox_cpuklt
        gadgetron_toolbox_mri_core
        )

target_include_directories(gadgetron_debugging
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

install(FILES
        gadgetron_debugging_export.h
        ${gadgetron_debugging_header_files}
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES ${gadgetron_debugging_config_files} DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

install(TARGETS gadgetron_debugging
        EXPORT gadgetron-export
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        COMPONENT main
)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)
