
set(EXAMPLE_SOURCES
        ImageLayerer.cpp
        ImageLayerer.h
        ImageInverter.cpp
        ImageInverter.h AcquisitionWaveformBranch.cpp AcquisitionWaveformBranch.h)

add_library(gadgetron_examples SHARED ${EXAMPLE_SOURCES})
set_target_properties(gadgetron_examples PROPERTIES
        VERSION ${GADGETRON_VERSION_STRING}
        SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_examples
        gadgetron_core
        gadgetron_core_parallel
        gadgetron_mricore
        gadgetron_toolbox_cpucore)

install(TARGETS gadgetron_examples DESTINATION lib COMPONENT main)

install(FILES
        config/external_example.xml
        config/parallel_bypass_example.xml
        DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
