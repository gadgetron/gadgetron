if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_FATWATER__)
endif (WIN32)

include_directories(
  ${Boost_INCLUDE_DIR} 
  ${ARMADILLO_INCLUDE_DIRS} 
  ${MKL_INCLUDE_DIR} 
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
  )

add_library(gadgetron_toolbox_fatwater SHARED 
  fatwater_export.h 
  fatwater.h
  fatwater.cpp
  )

set_target_properties(gadgetron_toolbox_fatwater PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_fatwater
  gadgetron_toolbox_log
  gadgetron_toolbox_cpucore 
  gadgetron_toolbox_cpucore_math 
  ${ARMADILLO_LIBRARIES} 
  )

install(TARGETS gadgetron_toolbox_fatwater DESTINATION lib COMPONENT main)

install(FILES
  fatwater_export.h 
  fatwater.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
