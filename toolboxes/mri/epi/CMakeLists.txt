    if (WIN32)
        ADD_DEFINITIONS(-D__BUILD_GADGETRON_EPI__)
        link_directories(${Boost_LIBRARY_DIRS})
    endif (WIN32)

    add_library(gadgetron_toolbox_epi ${LIBTYPE}
     EPIExport.h
     EPIReconXObject.h
     EPIReconXObjectFlat.h
     EPIReconXObjectTrapezoid.h)

    set_target_properties(gadgetron_toolbox_epi PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
    set_target_properties(gadgetron_toolbox_epi PROPERTIES LINKER_LANGUAGE CXX)

    target_link_libraries(gadgetron_toolbox_epi gadgetron_toolbox_cpucore gadgetron_toolbox_cpucore_math ${ARMADILLO_LIBRARIES})

    install (FILES 
        EPIExport.h 
        EPIReconXObject.h
        EPIReconXObjectFlat.h
        EPIReconXObjectTrapezoid.h
        DESTINATION include COMPONENT main)

    # install(TARGETS epi DESTINATION lib)

    #if (GTEST_FOUND)
    #    add_subdirectory(ut)
    #endif (GTEST_FOUND)

