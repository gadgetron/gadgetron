add_library(pingvin_gpuparallelmri SHARED
    GenericReconJob.h
    gpuCgKtSenseGadget.h
    gpuCgSenseGadget.h
    gpuCgSpiritGadget.h
    gpuGenericSensePrepGadget.h
    gpuSbSenseGadget.h
    gpuCgSenseGadget.cpp
    gpuCgKtSenseGadget.cpp
    gpuSbSenseGadget.cpp
    gpuGenericSensePrepGadget.cpp
    gpuCgSpiritGadget.cpp
    gpuSenseGadget.cpp
    gpuBufferSensePrepGadget.cpp
    gpuOsSenseGadget.cpp
    gpuNlcgSenseGadget.cpp
    gpuLALMSenseGadget.cpp
  )

set_target_properties(pingvin_gpuparallelmri PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(pingvin_gpuparallelmri
  pingvin_core
  pingvin_mricore
  pingvin_toolbox_log
  pingvin_toolbox_cpucore
  pingvin_toolbox_gpu)

target_link_libraries(pingvin_gpuparallelmri pingvin_toolbox_gpu)

target_include_directories(pingvin_gpuparallelmri
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})



install(TARGETS pingvin_gpuparallelmri
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install (FILES  GenericReconJob.h
                gpuCgKtSenseGadget.h
                gpuCgSenseGadget.h
                gpuCgSpiritGadget.h
                gpuGenericSensePrepGadget.h
                gpuSbSenseGadget.h
                gpuSenseGadget.h
                gpuBufferSensePrepGadget.h
                gpuOsSenseGadget.h
                gpuLALMSenseGadget.h
                gpuNlcgSenseGadget.h
                DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)

add_subdirectory(config)