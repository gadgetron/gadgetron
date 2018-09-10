if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_FATWATER__)
endif ()




include_directories(
  ${Boost_INCLUDE_DIR} 
  ${ARMADILLO_INCLUDE_DIRS}
  ${FFTW3_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_BINARY_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
  )

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
  ${ARMADILLO_LIBRARIES}
  ${Boost_LIBRARIES}
  )

install(TARGETS gadgetron_toolbox_fatwater DESTINATION lib COMPONENT main)

install(FILES
  fatwater_export.h 
  fatwater.h 
  correct_frequency_shift.h 
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
