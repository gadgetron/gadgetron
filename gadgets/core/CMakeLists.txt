if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_SPIRAL__)
endif ()

add_library(gadgetron_core SHARED
	Channel.h
	Channel.hpp
    Channel.cpp
    Context.cpp
    Stream.cpp
    MergeNode.cpp
    BranchNode.cpp)

target_link_libraries(gadgetron_core
		gadgetron_toolbox_cpucore)

target_include_directories(
		gadgetron_core
		PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}
    ${ISMRMRD_INCLUDE_DIR}
)
install(FILES
	Channel.h
        Channel.hpp
        Message.h
        Gadget.h
        Context.h
        Gadget.h
        Stream.h
        Node.h
        MergeNode.h
        BranchNode.h
		PropteryMixin.h
		Channel.h
		gadgetcore_export.h
		GadgetronExport.h
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES
        config/Example_Branching_Chain.xml
        DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
