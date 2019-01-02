configure_file(gadgetron_config.in gadgetron_config.h)

if (WIN32)
    add_definitions(-DCROW_MSVC_WORKAROUND)
endif ()


add_executable(gadgetron
        main.cpp
        GadgetronExport.h
        gadgetron_xml.h
        gadgetron_home.h
        gadgetron_home.cpp
        connection/Config.cpp
        connection/Config.h Server.cpp Server.h Connection.cpp Connection.h paths.cpp paths.h gadgetron_system_info.cpp connection/ConfigConnection.cpp connection/ConfigConnection.h connection/StreamConnection.cpp connection/StreamConnection.h connection/Handlers.cpp connection/Handlers.h connection/Writers.cpp connection/Writers.h connection/VoidConnection.cpp connection/VoidConnection.h connection/HeaderConnection.cpp connection/HeaderConnection.h connection/Loader.cpp connection/Loader.h)

target_link_libraries(gadgetron
        gadgetron_core
        gadgetron_toolbox_log
        pugixml
        boost
        ${ISMRMRD_LIBRARIES})

target_include_directories(gadgetron
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}
        )


if (BUILD_PYTHON_SUPPORT)
    add_definitions("-DCOMPILING_WITH_PYTHON_SUPPORT")
endif ()


if (CUDA_FOUND)
    target_link_libraries(gadgetron ${CUDA_LIBRARIES})
endif ()


if (GPERFTOOLS_PROFILER)
    message("Adding gperftools cpu profiler to gadgetron link assemblage.")
    target_link_libraries(gadgetron ${GPERFTOOLS_PROFILER} ${GPERFTOOLS_TCMALLOC})
endif ()


install(TARGETS gadgetron DESTINATION bin COMPONENT main)

