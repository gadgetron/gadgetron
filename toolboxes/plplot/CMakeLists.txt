# gadgetron_toolbox_plplot is a toolbox for plotting gadgetron data structures using PLplot

if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_TOOLBOX_PLPLOT__)
endif ()

include_directories(${PLPLOT_INCLUDE_DIR})
link_directories(${PLPLOT_LIB_DIR})


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

target_include_directories(gadgetron_toolbox_plplot
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS gadgetron_toolbox_plplot
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install (FILES  
        PLplotExport.h 
        ${header_files} 
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

if (GTEST_FOUND)
    add_subdirectory(ut)
endif ()

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)