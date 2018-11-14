
set(core_header_files
    Channel.h
    Message.h
    Gadget.h)


install(FILES ${core_header_files} DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH}/core COMPONENT main)
