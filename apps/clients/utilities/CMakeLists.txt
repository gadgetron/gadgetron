set(Boost_NO_BOOST_CMAKE ON)

if(WIN32)
  find_package(Boost COMPONENTS thread system date_time chrono REQUIRED)
else()
  find_package(Boost COMPONENTS thread system REQUIRED)
endif()

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif()


add_executable(gt_alive gt_alive.cpp)
add_executable(gtdependencyquery gt_query.cpp DependencyQueryReader.h)

target_link_libraries(gt_alive gadgetron_toolbox_cpucore 
                               gadgetron_toolbox_gadgettools 
			       gadgetron_toolbox_log)

target_link_libraries(gtdependencyquery gadgetron_toolbox_cpucore 
                                        gadgetron_toolbox_gadgettools 
                    					gadgetron_toolbox_log
                                        boost
                                         )

install(TARGETS gt_alive gtdependencyquery DESTINATION bin COMPONENT main)
install(FILES DependencyQueryReader.h DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
