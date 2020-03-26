
set(image_header_fiels 
        image_export.h
        morphology.h )

set(image_src_fiels morphology.cpp )
  
add_library(gadgetron_toolbox_cpu_image SHARED 
    ${image_header_fiels} ${image_src_fiels} )

set_target_properties(gadgetron_toolbox_cpu_image PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_cpu_image
    gadgetron_toolbox_log
    gadgetron_toolbox_cpucore 
    gadgetron_toolbox_cpucore_math 
    )

target_include_directories(gadgetron_toolbox_cpu_image
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )

install(TARGETS gadgetron_toolbox_cpu_image
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES
  ${image_header_fiels} 
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)