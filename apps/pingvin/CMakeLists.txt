find_package(PugiXML REQUIRED)

configure_file(pingvin_config.in pingvin_config.h)

add_executable(pingvin
        main.cpp
        initialization.cpp
        initialization.h
        system_info.cpp

        Loader.cpp
        Loader.h

        Config.cpp
        Config.h

        ErrorHandler.h

        nodes/Processable.h
        nodes/Processable.cpp

        nodes/Stream.cpp
        nodes/Stream.h
        nodes/Parallel.cpp
        nodes/Parallel.h
        nodes/ParallelProcess.cpp
        nodes/ParallelProcess.h
        nodes/PureStream.cpp
        nodes/PureStream.h
        )

target_link_libraries(pingvin
        pingvin_core
        pingvin_toolbox_log
        pingvin_toolbox_mri_core
        Boost::system
        Boost::filesystem
        Boost::program_options
        ${PUGIXML_LIBRARIES}
        GTBLAS
        ${CMAKE_DL_LIBS})

target_include_directories(pingvin
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR})


if (REQUIRE_SIGNED_CONFIG)
    target_link_libraries(pingvin GTBabylon)
endif()

if (CUDA_FOUND)
    target_link_libraries(pingvin ${CUDA_LIBRARIES})
endif ()

if (GPERFTOOLS_PROFILER)
    message("Adding gperftools cpu profiler to Pingvin link assemblage.")
    target_link_libraries(pingvin ${GPERFTOOLS_PROFILER} ${GPERFTOOLS_TCMALLOC})
endif ()

install(TARGETS pingvin DESTINATION bin COMPONENT main)
