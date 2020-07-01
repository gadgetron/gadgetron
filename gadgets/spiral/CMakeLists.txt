

if (CUDA_FOUND)
set(SPIRAL_HEADER_FILES gpuSpiralSensePrepGadget.h gpuSpiralDeblurGadget.h)
set(SPIRAL_SRC_FILES  gpuSpiralSensePrepGadget.cpp gpuSpiralDeblurGadget.cpp)

add_library(gadgetron_spiral SHARED ${SPIRAL_HEADER_FILES} ${SPIRAL_SRC_FILES})

set_target_properties(gadgetron_spiral PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_spiral
  gadgetron_core
  gadgetron_mricore
  gadgetron_toolbox_mri_core
  gadgetron_toolbox_log
	gadgetron_toolbox_spiral
  gadgetron_toolbox_cpucore)

if (CUDA_FOUND)
    target_link_libraries(gadgetron_spiral
                                gadgetron_gpuparallelmri
                                gadgetron_toolbox_gpucore
                                gadgetron_toolbox_gpunfft 
                                gadgetron_toolbox_gpusolvers 
                                gadgetron_toolbox_gpuoperators
                                gadgetron_toolbox_gpuparallelmri
                                gadgetron_toolbox_deblur )
endif()

install(TARGETS gadgetron_spiral
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install (FILES ${SPIRAL_HEADER_FILES}
                     DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

add_subdirectory(config)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)
endif()