

add_library(gadgetron_grappa_new SHARED
        ImageAccumulator.cpp
        WeightsCalculator.cpp
        Reconstruction.cpp
        AcquisitionFanout.cpp
        SliceAccumulator.cpp
        common/AcquisitionBuffer.cpp)

set_target_properties(gadgetron_grappa_new PROPERTIES
        VERSION ${GADGETRON_VERSION_STRING}
        SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_grappa_new
        gadgetron_core
        gadgetron_core_parallel
        gadgetron_toolbox_cpufft)

install(TARGETS gadgetron_grappa_new DESTINATION lib COMPONENT main)

install(FILES
        config/grappa.xml
        config/grappa_cpu.xml
        config/grappa_float.xml
        config/grappa_float_cpu.xml
        config/grappa_unoptimized.xml
        config/grappa_unoptimized_float.xml
        DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

