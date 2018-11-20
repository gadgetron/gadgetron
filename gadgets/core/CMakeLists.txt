include_directories(
    ${ISMRMRD_INCLUDE_DIR}
)

add_library(gadgetron_core SHARED
    Channel.cpp
    Context.cpp
    Stream.cpp
    MergeNode.cpp
    BranchNode.cpp)

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
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES
        config/Example_Branching_Chain.xml
        DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
