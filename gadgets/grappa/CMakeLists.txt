if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_GRAPPA__)
endif ()


add_library(gadgetron_grappa SHARED
        gadgetron_grappa_export.h
#        GrappaCalibrationBuffer.h
#        GrappaGadget.h
        GrappaUnmixingGadget.h
#        GrappaWeights.h
#        GrappaWeightsCalculator.h
#        GrappaGadget.cpp
#        GrappaCalibrationBuffer.cpp
#        GrappaWeights.cpp
#        GrappaWeightsCalculator.cpp
        GrappaUnmixingGadget.cpp
        )

set_target_properties(gadgetron_grappa PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
target_link_libraries(gadgetron_grappa
        gadgetron_core
        gadgetron_toolbox_log
        gadgetron_toolbox_gpuparallelmri
        gadgetron_toolbox_cpucore
        gadgetron_toolbox_cpucore_math
        gadgetron_toolbox_cpufft
        gadgetron_mricore
        )
if (CUDA_FOUND)
    target_link_libraries(gadgetron_grappa
            gadgetron_toolbox_gpufft
            )

endif ()

install(FILES gadgetron_grappa_export.h
        GrappaCalibrationBuffer.h
        GrappaGadget.h
        GrappaUnmixingGadget.h
        GrappaWeights.h
        GrappaWeightsCalculator.h
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_grappa DESTINATION lib COMPONENT main)

add_subdirectory(config)
