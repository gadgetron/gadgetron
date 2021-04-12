
add_library(gadgetron_session
        SessionServer.cpp
        SessionServer.h
        )
link_libraries(gadgetron_session
        Boost::beast
        )