if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_METRICS__)
endif ()

find_package(PrometheusCpp)
if (PROMETHEUSCPP_FOUND)
	message("Adding prometheus definition")
	add_definitions(-DPROMETHEUS_METRICS)
else()
	message("Prometheus metrics definition not added")
endif()

add_library(gadgetron_toolbox_metrics SHARED metrics.cpp)
target_include_directories(gadgetron_toolbox_metrics
		PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
set_target_properties(gadgetron_toolbox_metrics PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_metrics ${PROMETHEUSCPP_LIBRARIES})

install(TARGETS gadgetron_toolbox_metrics
		EXPORT gadgetron-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES metrics.h metrics_export.h DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)