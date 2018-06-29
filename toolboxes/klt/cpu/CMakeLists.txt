if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_CPUKLT__)
endif ()

include_directories(
  ${Boost_INCLUDE_DIR} 
  ${ARMADILLO_INCLUDE_DIRS} 
  ${MKL_INCLUDE_DIR} 
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_BINARY_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
  )

add_library(gadgetron_toolbox_cpuklt SHARED 
  cpuklt_export.h 
  hoNDKLT.h
  hoNDKLT.cpp
  )

set_target_properties(gadgetron_toolbox_cpuklt PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_cpuklt
                    gadgetron_toolbox_cpucore
                    gadgetron_toolbox_log
                    gadgetron_toolbox_cpucore_math
                    )

install(TARGETS gadgetron_toolbox_cpuklt DESTINATION lib COMPONENT main)

install(FILES
  cpuklt_export.h 
  hoNDKLT.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
