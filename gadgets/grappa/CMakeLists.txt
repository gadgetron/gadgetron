set(RTGRAPPA_SOURCES
        ImageAccumulator.cpp
        WeightsCalculator.cpp
        Unmixing.cpp
        AcquisitionFanout.cpp
        SliceAccumulator.cpp
        common/AcquisitionBuffer.cpp
        common/AcquisitionBuffer.h
        common/grappa_common.h
        common/AnnotatedAcquisition.h
        ChannelReorderer.cpp
        ChannelReorderer.h
        cpu/WeightsCore.cpp
        cpu/WeightsCore.h)

if (CUDA_FOUND)
    set(RTGRAPPA_SOURCES ${RTGRAPPA_SOURCES}
            gpu/WeightsCore.cpp
            gpu/WeightsCore.h)
endif()

add_library(pingvin_grappa SHARED ${RTGRAPPA_SOURCES})
set_target_properties(pingvin_grappa PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(pingvin_grappa
        pingvin_core
        pingvin_core_parallel
        pingvin_mricore
        pingvin_toolbox_cpucore
        pingvin_toolbox_cpucore_math
        pingvin_toolbox_cpufft)

if (CUDA_FOUND)
    target_link_libraries(pingvin_grappa
            pingvin_toolbox_gpu
    )
endif()

install(TARGETS pingvin_grappa
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
        config/grappa.xml
        config/grappa_cpu.xml
        config/grappa_float.xml
        config/grappa_float_cpu.xml
        config/grappa_unoptimized.xml
        config/grappa_unoptimized_float.xml
        DESTINATION ${PINGVIN_INSTALL_CONFIG_PATH} COMPONENT main)