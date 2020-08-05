if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_FATWATER__)
endif ()

add_library(gadgetron_toolbox_fatwater SHARED
  fatwater_export.h 
  fatwater.h
  fatwater.cpp
        graph_cut.cpp ImageGraph.cpp correct_frequency_shift.h correct_frequency_shift.cpp bounded_field_map.cpp)

set_target_properties(gadgetron_toolbox_fatwater PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_fatwater
  gadgetron_toolbox_log
  gadgetron_toolbox_cpucore 
  gadgetron_toolbox_cpucore_math
  gadgetron_toolbox_cpufft
  armadillo
  )

target_include_directories(
    gadgetron_toolbox_fatwater
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
  )

install(TARGETS gadgetron_toolbox_fatwater
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
  fatwater_export.h 
  fatwater.h 
  correct_frequency_shift.h 
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)