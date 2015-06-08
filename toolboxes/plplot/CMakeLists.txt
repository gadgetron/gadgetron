# gadgetron_toolbox_plplot is a toolbox for plotting gadgetron data structures using PLplot

if (WIN32)
    ADD_DEFINITIONS(-D__BUILD_GADGETRON_TOOLBOX_PLPLOT__)
endif (WIN32)

include_directories(${PLPLOT_INCLUDE_DIR})
link_directories(${PLPLOT_LIB_DIR})

include_directories(
    ${Boost_INCLUDE_DIR} 
    ${ARMADILLO_INCLUDE_DIRS} 
    ${ISMRMRD_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
)

set(header_files GtPLplot.h)
set(src_files GtPLplot.cpp)

set(plplot_files PLplotExport.h 
            ${header_files} 
            ${src_files} )

source_group(plplot FILES ${header_files} ${src_files})

add_library(gadgetron_toolbox_plplot SHARED 
            ${plplot_files} )

set_target_properties(gadgetron_toolbox_plplot PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_plplot 
                      ${PLPLOT_LIBRARIES}
                      gadgetron_toolbox_cpucore
                      gadgetron_toolbox_log
                      gadgetron_toolbox_cpucore_math 
                      gadgetron_toolbox_mri_core
                      gadgetron_toolbox_cpucore_math )

install(TARGETS gadgetron_toolbox_plplot DESTINATION lib COMPONENT main)

install (FILES  
        PLplotExport.h 
        ${header_files} 
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
