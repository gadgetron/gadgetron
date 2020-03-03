


add_library(gadgetron_toolbox_t1 SHARED
        t1fit.cpp)

set_target_properties(gadgetron_toolbox_t1 PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_t1
        gadgetron_toolbox_log
        gadgetron_toolbox_cpucore
        gadgetron_toolbox_cpucore_math
		gadgetron_toolbox_cpu_solver
        armadillo
        )
target_include_directories(gadgetron_toolbox_t1
        PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        )

install(TARGETS gadgetron_toolbox_t1
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES
        t1fit.h
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
