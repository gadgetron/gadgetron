set(Boost_NO_BOOST_CMAKE ON)

if(WIN32)
  find_package(Boost COMPONENTS thread system date_time chrono REQUIRED)
else(WIN32)
  find_package(Boost COMPONENTS thread system REQUIRED)
endif(WIN32)

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories(
    ${CMAKE_SOURCE_DIR}/apps/gadgetron
    ${CMAKE_BINARY_DIR}/apps/gadgetron
    ${CMAKE_SOURCE_DIR}/gadgets/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
    ${ACE_INCLUDE_DIR} 
    ${Boost_INCLUDE_DIR} 
    ${ISMRMRD_INCLUDE_DIR}
    )

add_executable(gt_alive gt_alive.cpp)
add_executable(gtdependencyquery gt_query.cpp DependencyQueryReader.h)

target_link_libraries(gt_alive gadgetron_toolbox_cpucore 
                               gadgetron_toolbox_gadgettools 
			       gadgetron_toolbox_log
                               optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
                               ${Boost_LIBRARIES} 
                               ${ISMRMRD_LIBRARIES} )

target_link_libraries(gtdependencyquery gadgetron_toolbox_cpucore 
                                        gadgetron_toolbox_gadgettools 
					gadgetron_toolbox_log
                                        optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
                                        ${Boost_LIBRARIES} 
                                        ${ISMRMRD_LIBRARIES} )

install(TARGETS gt_alive gtdependencyquery DESTINATION bin COMPONENT main)
install(FILES DependencyQueryReader.h DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
