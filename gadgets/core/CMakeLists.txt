
add_library(gadgetron_core SHARED
    Gadget.cpp
    Node.cpp
    ConvergeNode.cpp
    DisperseNode.cpp)

install(FILES
        Gadget.h
        Node.h
        ConvergeNode.h
        DisperseNode.h
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES
        config/Example_Branching_Chain.xml
        DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
