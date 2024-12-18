
add_library(pingvin_toolbox_spiral SHARED TrajectoryParameters.cpp vds.cpp TrajectoryParameters.h vds.h)


target_link_libraries(pingvin_toolbox_spiral pingvin_toolbox_cpucore pingvin_toolbox_mri_core)

set_target_properties(pingvin_toolbox_spiral PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_include_directories(
        pingvin_toolbox_spiral
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

install(TARGETS pingvin_toolbox_spiral
        EXPORT pingvin-export
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        COMPONENT main
        )

install(FILES
        vds.h
        TrajectoryParameters.h
        DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)