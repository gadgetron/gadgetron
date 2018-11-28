if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_MRI_IMAGE__)
endif ()


set( mri_image_header_files
        mri_image_export.h
        hoMRImage.h 
        hoMRImage.hxx )

set( mri_image_source_files )

add_library(gadgetron_toolbox_mri_image SHARED ${mri_image_header_files} ${mri_image_source_files} )

set_target_properties(gadgetron_toolbox_mri_image PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
set_target_properties(gadgetron_toolbox_mri_image PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(gadgetron_toolbox_mri_image
                    gadgetron_toolbox_cpucore
                    gadgetron_toolbox_cpucore_math
                     )
target_include_directories(gadgetron_toolbox_mri_image
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})

# install(TARGETS gadgetron_toolbox_mri_image DESTINATION lib COMPONENT main)

install(FILES ${mri_image_header_files}
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
