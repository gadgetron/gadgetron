include_directories(${Boost_INCLUDE_DIR})

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math)

add_library(gadgetron_toolbox_mri_core SHARED
  mri_core.cpp)

target_link_libraries(gadgetron_toolbox_mri_core
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_cpucore_math)

set_target_properties(gadgetron_toolbox_mri_core PROPERTIES COMPILE_DEFINITIONS "__BUILD_GADGETRON_MRI_CORE__")
set_target_properties(gadgetron_toolbox_mri_core PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

install(TARGETS gadgetron_toolbox_mri_core DESTINATION lib COMPONENT main)

install(FILES 
  mri_core.h
  mri_core_data.h
  mri_core_export.h
  DESTINATION include COMPONENT main)
