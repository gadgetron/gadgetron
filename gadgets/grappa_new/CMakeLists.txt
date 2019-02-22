

add_library(gadgetron_grappa_new SHARED
        ImageAccumulator.cpp
        WeightsCalculator.cpp
        Unmixing.cpp
        AcquisitionFanout.cpp
        SliceAccumulator.cpp
        common/AcquisitionBuffer.cpp common/grappa_common.h ChannelReorderer.cpp ChannelReorderer.h common/AnnotatedAcquisition.h CPU/WeightsCore.cpp CPU/WeightsCore.h GPU/WeightsCore.cpp GPU/WeightsCore.h)

set_target_properties(gadgetron_grappa_new PROPERTIES
        VERSION ${GADGETRON_VERSION_STRING}
        SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_grappa_new
        gadgetron_core
        gadgetron_core_parallel
        gadgetron_mricore
        gadgetron_toolbox_cpucore
        gadgetron_toolbox_cpucore_math
        gadgetron_toolbox_cpufft)

install(TARGETS gadgetron_grappa_new DESTINATION lib COMPONENT main)

install(FILES
#        config/grappa.xml
        config/grappa_cpu.xml
#        config/grappa_float.xml
        config/grappa_float_cpu.xml
#        config/grappa_unoptimized.xml
#        config/grappa_unoptimized_float.xml
        config/grappa_uncombined_channel_example.xml
        DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

