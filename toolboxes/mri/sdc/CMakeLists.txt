add_subdirectory(cpu)


add_library(pingvin_toolbox_sdc INTERFACE)
target_include_directories(pingvin_toolbox_sdc
        INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR})

install(FILES
    SDC.h
    DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)