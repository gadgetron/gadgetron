    if (WIN32)
        add_definitions(-D__BUILD_GADGETRON_EPI__)
    endif ()

    add_library(gadgetron_toolbox_epi SHARED
            EPIExport.h
            EPIReconXObject.h
            EPIReconXObjectFlat.h
            EPIReconXObjectTrapezoid.h)

    set_target_properties(gadgetron_toolbox_epi PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
    set_target_properties(gadgetron_toolbox_epi PROPERTIES LINKER_LANGUAGE CXX)

    target_link_libraries(gadgetron_toolbox_epi gadgetron_toolbox_cpucore gadgetron_toolbox_cpucore_math gadgetron_toolbox_log boost)

    target_include_directories(gadgetron_toolbox_epi
            PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR})

    install(FILES
            EPIExport.h
            EPIReconXObject.h
            EPIReconXObjectFlat.h
            EPIReconXObjectTrapezoid.h
            DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

    # install(TARGETS epi DESTINATION lib)

    #if (GTEST_FOUND)
    #    add_subdirectory(ut)
    #endif ()

