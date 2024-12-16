add_library(pingvin_gpuradial SHARED
  gpuRadialPrepGadget.h gpuRadialPrepGadget.cpp
  gpuRadialSensePrepGadget.h gpuRadialSensePrepGadget.cpp
  gpuRadialSpiritPrepGadget.h gpuRadialSpiritPrepGadget.cpp
  gpuRetroGatedSensePrepGadget.h gpuRetroGatedSensePrepGadget.cpp
  RadialPhaseCorrectionGadget.h RadialPhaseCorrectionGadget.cpp
  )

set_target_properties(pingvin_gpuradial PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(pingvin_gpuradial
        pingvin_core
        pingvin_mricore
        pingvin_toolbox_log
        pingvin_gpuparallelmri
        pingvin_toolbox_cpucore
        pingvin_toolbox_cpucore_math
        pingvin_toolbox_gpu)

install (FILES  gpuRadialPrepGadget.h
                gpuRadialSensePrepGadget.h
                gpuRadialSpiritPrepGadget.h
                gpuRetroGatedSensePrepGadget.h
                RadialPhaseCorrectionGadget.h
                DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS pingvin_gpuradial
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

add_subdirectory(config)