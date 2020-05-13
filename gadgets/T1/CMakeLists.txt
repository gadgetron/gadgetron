add_library(gadgetron_t1 SHARED
        T1MocoGadget.cpp)

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
