if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_METRICS__)
endif ()

find_package(PrometheusCpp)
if (PROMETHEUSCPP_FOUND)
	add_definitions(-DPROMETHEUS_METRICS)
	add_library(gadgetron_toolbox_metrics SHARED Metrics.cpp PrometheusMetricsSink.cpp)
else()
	add_library(gadgetron_toolbox_metrics SHARED Metrics.cpp)
endif()

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

set(metrics_header_files
	Metrics.h
	metrics_export.h
	MetricsSink.h
	MetricsSinkFactory.h
	NoOpMetricsSink.h
	PrometheusMetricsSink.h
)

install(FILES ${metrics_header_files} DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)