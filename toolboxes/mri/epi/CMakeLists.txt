    add_library(pingvin_toolbox_epi  INTERFACE)

    target_link_libraries(pingvin_toolbox_epi INTERFACE pingvin_toolbox_cpucore pingvin_toolbox_cpucore_math pingvin_toolbox_log )

    target_include_directories(pingvin_toolbox_epi
            INTERFACE
            ${CMAKE_CURRENT_SOURCE_DIR})

    install(FILES
            EPIReconXObject.h
            EPIReconXObjectFlat.h
            EPIReconXObjectTrapezoid.h
            DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)

    # install(TARGETS epi DESTINATION lib)

    #if (GTEST_FOUND)
    #    add_subdirectory(ut)
    #endif ()