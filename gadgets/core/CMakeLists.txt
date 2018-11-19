include_directories(
    ${ISMRMRD_INCLUDE_DIR}
)
set(core_header_files
    Channel.h
    Channel.hpp
    Message.h
    Gadget.h)

set(core_source_files
    Channel.cpp)

install(FILES ${core_header_files} DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH}/core COMPONENT main)

add_library(core ${core_header_files} ${core_source_files} )
