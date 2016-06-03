IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETS_FATWATER__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
  ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
  ${CMAKE_SOURCE_DIR}/toolboxes/fatwater
)

add_library(gadgetron_fatwater SHARED 
  FatWaterGadget.h FatWaterGadget.cpp
)

set_target_properties(gadgetron_fatwater PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(
  gadgetron_fatwater
  gadgetron_gadgetbase
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_cpufft
  gadgetron_toolbox_cpucore_math
  gadgetron_toolbox_log
  gadgetron_toolbox_fatwater
  ${ISMRMRD_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
)

install(FILES 
  FatWaterGadget.h
  gadgetron_fatwater_export.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_fatwater DESTINATION lib COMPONENT main)

install(FILES
  config/Generic_Cartesian_Grappa_FatWater_Proto.xml
  DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
