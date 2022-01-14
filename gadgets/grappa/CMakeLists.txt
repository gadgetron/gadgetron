

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

add_library(gadgetron_grappa SHARED ${RTGRAPPA_SOURCES})
#set_target_properties(gadgetron_grappa PROPERTIES
#        VERSION ${GADGETRON_VERSION_STRING}
#        SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_grappa
        gadgetron_core
        gadgetron_core_parallel
        gadgetron_mricore
        gadgetron_toolbox_cpucore
        gadgetron_toolbox_cpucore_math
        gadgetron_toolbox_cpufft)

if (CUDA_FOUND)
    target_link_libraries(gadgetron_grappa
            gadgetron_toolbox_gpu
    )
endif()

install(TARGETS gadgetron_grappa
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
        DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)