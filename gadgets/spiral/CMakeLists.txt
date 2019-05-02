if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_GADGETCORE__)
endif ()

set(SPIRAL_HEADER_FILES gadgetron_spiral_export.h SpiralToGenericGadget.h vds.h TrajectoryParameters.h)
set(SPIRAL_SRC_FILES vds.cpp  SpiralToGenericGadget.cpp TrajectoryParameters.cpp)

if (CUDA_FOUND)
set(SPIRAL_HEADER_FILES ${SPIRAL_HEADER_FILES} gpuSpiralSensePrepGadget.h)
set(SPIRAL_SRC_FILES ${SPIRAL_SRC_FILES}  gpuSpiralSensePrepGadget.cpp)
endif()

add_library(gadgetron_spiral SHARED ${SPIRAL_HEADER_FILES} ${SPIRAL_SRC_FILES})

set_target_properties(gadgetron_spiral PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_spiral
  gadgetron_core
  gadgetron_mricore
  gadgetron_toolbox_mri_core
  gadgetron_toolbox_log
  gadgetron_toolbox_cpucore)

if (CUDA_FOUND)
    target_link_libraries(gadgetron_spiral
                                gadgetron_gpuparallelmri
                                gadgetron_toolbox_gpucore
                                gadgetron_toolbox_gpunfft 
                                gadgetron_toolbox_gpusolvers 
                                gadgetron_toolbox_gpuoperators
                                gadgetron_toolbox_gpuparallelmri )
endif()


install (TARGETS gadgetron_spiral DESTINATION lib COMPONENT main)
install (FILES ${SPIRAL_HEADER_FILES}
                     DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

add_subdirectory(config)
