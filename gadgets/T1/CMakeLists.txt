
set(gadgetron_t1_config_files config/MOLLI_T1_Moco.xml config/MOLLI_T1_Moco_istore.xml)

source_group(config FILES ${gadgetron_t1_config_files})

add_library(gadgetron_t1 SHARED
        T1MocoGadget.cpp T1_utils.hpp)

set_target_properties(gadgetron_t1 PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_t1
    gadgetron_core
    gadgetron_toolbox_log
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_cpufft
    gadgetron_toolbox_demons
    gadgetron_toolbox_t1
    gadgetron_toolbox_mri_core
)

install(TARGETS gadgetron_t1
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES ${gadgetron_t1_config_files} DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
