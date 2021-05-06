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
        connection/config/Config.cpp
        connection/config/Config.h
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
        connection/nodes/Stream.cpp
        connection/nodes/Stream.h
        connection/nodes/Parallel.cpp
        connection/nodes/Parallel.h
        connection/nodes/Distributed.cpp
        connection/nodes/Distributed.h
        connection/nodes/External.cpp
        connection/nodes/External.h
        connection/core/Processable.h
        connection/nodes/common/ExternalChannel.cpp
        connection/nodes/common/ExternalChannel.h
        connection/nodes/external/Matlab.cpp
        connection/nodes/external/Matlab.h
        connection/nodes/external/Python.cpp
        connection/nodes/external/Python.h
        connection/nodes/common/Discovery.cpp
        connection/nodes/ParallelProcess.cpp
        connection/nodes/ParallelProcess.h
        connection/nodes/PureStream.cpp
        connection/nodes/PureStream.h
        connection/nodes/PureDistributed.cpp
        connection/nodes/PureDistributed.h
        connection/nodes/common/External.h
        connection/nodes/common/External.cpp
        connection/nodes/common/Serialization.cpp
        connection/nodes/common/Serialization.h
        connection/nodes/common/Configuration.cpp
        connection/nodes/common/Configuration.h
        connection/nodes/distributed/Pool.h
        connection/nodes/distributed/Worker.cpp
        connection/nodes/distributed/Worker.h
        connection/nodes/common/Closer.h
        connection/nodes/distributed/Pool.cpp connection/core/Processable.cpp)

target_link_libraries(gadgetron
        gadgetron_core
        gadgetron_session
        gadgetron_toolbox_log
        Boost::system
        Boost::filesystem
        Boost::program_options
        GTBLAS
        ${CMAKE_DL_LIBS})

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



