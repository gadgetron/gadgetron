find_package(Boost COMPONENTS date_time filesystem system thread regex REQUIRED)
find_package( Threads )

if (WIN32)
  add_definitions(-DCROW_MSVC_WORKAROUND)
endif(WIN32)

include_directories(${Boost_INCLUDE_DIR})

add_library(gadgetron_toolbox_rest SHARED
  gadgetron_rest.cpp
  gadgetron_rest.h
  gadgetron_rest_exports.h
)

target_link_libraries(gadgetron_toolbox_rest
		      gadgetron_toolbox_log
		      ${Boost_LIBRARIES}
		      ${CMAKE_THREAD_LIBS_INIT})

set_target_properties(gadgetron_toolbox_rest PROPERTIES COMPILE_DEFINITIONS "__BUILD_GADGETRON_REST__")
set_target_properties(gadgetron_toolbox_rest  PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

install(TARGETS gadgetron_toolbox_rest DESTINATION lib COMPONENT main)

install(FILES 
  gadgetron_rest.h
  gadgetron_rest_exports.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
