

set(SPIRAL_HEADER_FILES SpiralToGenericGadget.h)
set(SPIRAL_SRC_FILES SpiralToGenericGadget.cpp)

if (CUDA_FOUND)
set(SPIRAL_HEADER_FILES gpuSpiralSensePrepGadget.h gpuSpiralDeblurGadget.h ${SPIRAL_HEADER_FILES})
set(SPIRAL_SRC_FILES  gpuSpiralSensePrepGadget.cpp gpuSpiralDeblurGadget.cpp ${SPIRAL_SRC_FILES})
endif()

add_library(pingvin_spiral SHARED ${SPIRAL_HEADER_FILES} ${SPIRAL_SRC_FILES})

set_target_properties(pingvin_spiral PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(pingvin_spiral
  pingvin_core
  pingvin_mricore
  pingvin_toolbox_mri_core
  pingvin_toolbox_log
	pingvin_toolbox_spiral
  pingvin_toolbox_cpucore)

if (CUDA_FOUND)
    target_link_libraries(pingvin_spiral
                                pingvin_gpuparallelmri
                                pingvin_toolbox_gpu
                                pingvin_core
                                pingvin_toolbox_deblur )
endif()

install(TARGETS pingvin_spiral
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install (FILES ${SPIRAL_HEADER_FILES}
                     DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)

add_subdirectory(config)