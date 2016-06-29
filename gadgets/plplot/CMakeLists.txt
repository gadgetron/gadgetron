if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_PLPLOT__)
endif ()

include_directories(${PLPLOT_INCLUDE_DIR})
link_directories(${PLPLOT_LIB_DIR})

include_directories(
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/log
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/plplot
    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
    ${ARMADILLO_INCLUDE_DIRS}
    ${MKL_INCLUDE_DIR}
    ${ISMRMRD_INCLUDE_DIR}
)

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
    gadgetron_gadgetbase
    gadgetron_mricore
    gadgetron_toolbox_log
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_plplot
    ${ISMRMRD_LIBRARIES} 
    optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
    ${BOOST_LIBRARIES}
)

install(FILES 
    gadgetron_plplot_export.h
    ${header_files}
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_plplot DESTINATION lib COMPONENT main)

install (FILES  ${config_files} 
        DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
