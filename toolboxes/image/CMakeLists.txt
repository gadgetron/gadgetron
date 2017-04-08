if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_IMAGE__)
endif ()

include_directories(
    ${Boost_INCLUDE_DIR} 
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/image
    )

set(image_header_fiels 
        image_export.h 
        morphology.h )

set(image_src_fiels morphology.cpp )
  
add_library(gadgetron_toolbox_image SHARED 
    ${image_header_fiels} ${image_src_fiels} )

set_target_properties(gadgetron_toolbox_image PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_image
    gadgetron_toolbox_log
    gadgetron_toolbox_cpucore 
    gadgetron_toolbox_cpucore_math 
    )

install(TARGETS gadgetron_toolbox_image DESTINATION lib COMPONENT main)

install(FILES
  ${image_header_fiels} 
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
