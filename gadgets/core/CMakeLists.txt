
add_library(gadgetron_core SHARED
    Context.cpp
    Gadget.cpp
    Stream.cpp
    Node.cpp
    MergeNode.cpp
    BranchNode.cpp)

install(FILES
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
