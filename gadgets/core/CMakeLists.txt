if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_CORE__)
endif ()

add_subdirectory(readers)
add_subdirectory(writers)
add_subdirectory(gadgets)

add_library(gadgetron_core SHARED
	Channel.h
	Channel.hpp
    Channel.cpp
    Context.cpp
	Gadget.cpp
    Reader.h
    Writer.h
	LegacyACE.cpp
	Message.cpp Writer.hpp)

target_link_libraries(gadgetron_core
		gadgetron_toolbox_cpucore)

target_include_directories(
    gadgetron_core
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${ISMRMRD_INCLUDE_DIR})


install(TARGETS gadgetron_core DESTINATION lib COMPONENT main)

install(FILES
    Channel.h
    Channel.hpp
    Message.h
    Gadget.h
    Context.h
    Gadget.h
    Reader.h
    Writer.h
    Node.h
    LegacyACE.h
    PropteryMixin.h
    Channel.h
    gadgetcore_export.h
    GadgetronExport.h
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES
        config/example_branching_chain.xml
        config/simple_chain.xml
    DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
