if (WIN32)
  add_definitions(-D__BUILD_GADGETS_FATWATER__)
endif ()


add_library(gadgetron_fatwater SHARED 
  FatWaterGadget.h FatWaterGadget.cpp
)

set_target_properties(gadgetron_fatwater PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(
  gadgetron_fatwater
  gadgetron_core
  gadgetron_mricore
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_cpufft
  gadgetron_toolbox_cpucore_math
  gadgetron_toolbox_log
  gadgetron_toolbox_fatwater
  ${ISMRMRD_LIBRARIES}
)

install(FILES 
  FatWaterGadget.h
  gadgetron_fatwater_export.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

  install(TARGETS gadgetron_fatwater
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES
        config/Generic_Cartesian_Grappa_FatWater.xml
  DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)