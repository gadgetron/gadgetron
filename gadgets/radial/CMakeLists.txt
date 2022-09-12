if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_RADIAL__)
endif ()

list(APPEND PHASE_GADGET RadialPhaseCorrectionGadget.h RadialPhaseCorrectionGadget.cpp)
install (FILES  RadialPhaseCorrectionGadget.h DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

add_library(gadgetron_gpuradial SHARED 
  gadgetron_radial_export.h
  gpuRadialPrepGadget.h gpuRadialPrepGadget.cpp 
  gpuRadialSensePrepGadget.h gpuRadialSensePrepGadget.cpp 
  gpuRadialSpiritPrepGadget.h gpuRadialSpiritPrepGadget.cpp 
  gpuRetroGatedSensePrepGadget.h gpuRetroGatedSensePrepGadget.cpp
  ${PHASE_GADGET})

set_target_properties(gadgetron_gpuradial PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_gpuradial
        gadgetron_core
        gadgetron_mricore
        gadgetron_toolbox_log
        gadgetron_gpuparallelmri
        gadgetron_toolbox_cpucore
        gadgetron_toolbox_cpucore_math
        gadgetron_toolbox_gpu)

install (FILES  gadgetron_radial_export.h
                gpuRadialPrepGadget.h 
                gpuRadialSensePrepGadget.h 
                gpuRadialSpiritPrepGadget.h 
                gpuRetroGatedSensePrepGadget.h 
                DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_gpuradial
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

add_subdirectory(config)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)