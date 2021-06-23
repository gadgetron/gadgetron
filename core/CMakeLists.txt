if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_CORE__)
endif ()

add_subdirectory(readers)
add_subdirectory(writers)
add_subdirectory(parallel)
add_subdirectory(distributed)

add_library(gadgetron_core SHARED
        Channel.cpp
        Gadget.cpp
        LegacyACE.cpp
        Message.cpp
        Response.cpp
        Storage.cpp
        Process.cpp
        gadgetron_paths.cpp
        io/from_string.cpp)
set_target_properties(gadgetron_core PROPERTIES
        VERSION ${GADGETRON_VERSION_STRING}
        SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_core
        ISMRMRD::ISMRMRD
        gadgetron_toolbox_cpucore
        Boost::boost
        Boost::filesystem
        )

target_include_directories(gadgetron_core PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
        )

install(TARGETS gadgetron_core
        EXPORT gadgetron-export
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
        Gadget.h
        Reader.h
        Types.h
        Types.hpp
        TypeTraits.h
        Writer.h
        Writer.hpp
        Node.h
        PureGadget.h
        LegacyACE.h
        PropertyMixin.h
        GadgetContainerMessage.h
        ChannelAlgorithms.h
        variant.hpp
        MessageID.h
        Storage.h
        Process.h
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES
        io/adapt_struct.h
        io/from_string.h
        io/ismrmrd_types.h
        io/primitives.h
        io/primitives.hpp
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH}/io COMPONENT main)
install(FILES
        config/distributed_default.xml
        config/distributed_generic_default.xml
        config/distributed_image_default.xml
        DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)