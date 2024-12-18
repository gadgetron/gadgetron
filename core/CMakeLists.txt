add_compile_options(-Wall -Werror)
add_subdirectory(parallel)

add_library(pingvin_core SHARED
        Channel.cpp
        Gadget.cpp
        Message.cpp
        Process.cpp
        io/from_string.cpp)

set_target_properties(pingvin_core PROPERTIES
        VERSION ${PINGVIN_VERSION_STRING}
        SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(pingvin_core
        mrd::mrd
        pingvin_toolbox_cpucore
        Boost::boost
        Boost::filesystem
        )

target_include_directories(pingvin_core PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
        )

install(TARGETS pingvin_core
        EXPORT pingvin-export
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        COMPONENT main
        )


install(FILES
        Channel.h
        Channel.hpp
        ChannelIterator.h
        Message.h
        Message.hpp
        MPMCChannel.h
        Gadget.h
        Context.h
        Node.h
        PureGadget.h
        PropertyMixin.h
        GadgetContainerMessage.h
        ChannelAlgorithms.h
        Process.h
        DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES
        io/from_string.h
        io/primitives.h
        io/primitives.hpp
        DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH}/io COMPONENT main)