if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_PLPLOT__)
endif ()

include_directories(${PLPLOT_INCLUDE_DIR})
link_directories(${PLPLOT_LIB_DIR})


set( header_files NoiseCovariancePlottingGadget.h )
set( src_files NoiseCovariancePlottingGadget.cpp )
set( config_files Generic_Cartesian_Grappa_SNR_CoilQA.xml )

source_group(config FILES ${config_files})

add_library(gadgetron_plplot SHARED 
    gadgetron_plplot_export.h 
    ${header_files} 
    ${src_files}
    ${config_files}
)

set_target_properties(gadgetron_plplot PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_plplot 
    gadgetron_core
    gadgetron_mricore
    gadgetron_toolbox_log
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_plplot
)

install(FILES 
    gadgetron_plplot_export.h
    ${header_files}
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_plplot
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install (FILES  ${config_files} 
        DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)