configure_file(gadgetron_config.in gadgetron_config.h)



add_executable(gadgetron
        main.cpp
        Server.cpp
        Server.h
        Connection.cpp
        Connection.h
        paths.cpp
        paths.h
        initialization.cpp
        initialization.h
        system_info.cpp
        connection/Config.cpp
        connection/Config.h
        connection/ConfigConnection.cpp
        connection/ConfigConnection.h
        connection/StreamConnection.cpp
        connection/StreamConnection.h
        connection/Handlers.cpp
        connection/Handlers.h
        connection/Writers.cpp
        connection/Writers.h
        connection/VoidConnection.cpp
        connection/VoidConnection.h
        connection/HeaderConnection.cpp
        connection/HeaderConnection.h
        connection/Loader.cpp
        connection/Loader.h
        connection/Core.cpp
        connection/Core.h
        connection/SocketStreamBuf.cpp
        connection/SocketStreamBuf.h
        connection/stream/Stream.cpp
        connection/stream/Stream.h
        connection/stream/Parallel.cpp
        connection/stream/Parallel.h
        connection/stream/Distributed.cpp
        connection/stream/Distributed.h
        connection/stream/External.cpp
        connection/stream/External.h
        connection/stream/Processable.h
        connection/stream/common/ExternalChannel.cpp
        connection/stream/common/ExternalChannel.h
        connection/stream/external/Matlab.cpp
        connection/stream/external/Matlab.h
        connection/stream/external/Python.cpp
        connection/stream/external/Python.h
        connection/stream/common/Discovery.cpp
        connection/stream/Processable.cpp
        connection/stream/ParallelProcess.cpp
        connection/stream/ParallelProcess.h
        connection/stream/PureStream.cpp
        connection/stream/PureStream.h
        connection/stream/PureDistributed.cpp
        connection/stream/PureDistributed.h
        connection/stream/common/External.h
        connection/stream/common/External.cpp
        connection/stream/common/Serialization.cpp
        connection/stream/common/Serialization.h
        connection/stream/common/Configuration.cpp
        connection/stream/common/Configuration.h
        connection/stream/distributed/Pool.h
        connection/stream/distributed/Worker.cpp
        connection/stream/distributed/Worker.h
        connection/stream/common/Closer.h
        connection/stream/distributed/Pool.cpp )

target_link_libraries(gadgetron
        gadgetron_core
        gadgetron_toolbox_log
        Boost::system
        Boost::filesystem
        Boost::program_options
        GTBLAS
        ${PUGIXML_LIBRARIES}
        ${CMAKE_DL_LIBS}
        )

target_include_directories(gadgetron
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR})


if (REQUIRE_SIGNED_CONFIG)
    target_link_libraries(gadgetron GTBabylon)
endif()

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



