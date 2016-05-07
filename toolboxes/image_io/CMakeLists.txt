
if (WIN32)
    ADD_DEFINITIONS(-D__BUILD_GADGETRON_IMAGE_ANALYZE_IO__)
endif (WIN32)

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories(
    ${ACE_INCLUDE_DIR} 
    ${Boost_INCLUDE_DIR}
    ${ISMRMRD_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/image_io
)

set( image_io_header_files 
     ImageIOExport.h 
     ImageIOBase.h
     ImageIOAnalyze.h)

set( image_io_src_files 
     ImageIOBase.cpp
     ImageIOAnalyze.cpp)

add_library(gadgetron_toolbox_image_analyze_io SHARED ${image_io_header_files} ${image_io_src_files} )
set_target_properties(gadgetron_toolbox_image_analyze_io PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
target_link_libraries(gadgetron_toolbox_image_analyze_io gadgetron_toolbox_log)

install(TARGETS gadgetron_toolbox_image_analyze_io DESTINATION lib COMPONENT main)

install (FILES  
        ${image_io_header_files} 
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
