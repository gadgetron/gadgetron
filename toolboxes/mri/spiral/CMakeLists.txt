
add_library(gadgetron_toolbox_spiral SHARED TrajectoryParameters.cpp vds.cpp TrajectoryParameters.h vds.h)


target_link_libraries(gadgetron_toolbox_spiral gadgetron_toolbox_cpucore gadgetron_toolbox_mri_core)

set_target_properties(gadgetron_toolbox_spiral PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_include_directories(
        gadgetron_toolbox_spiral
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

install(TARGETS gadgetron_toolbox_spiral
        EXPORT gadgetron-export
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        COMPONENT main
        )



install(FILES
        vds.h
        TrajectoryParameters.h
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)
