if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_MRI_IMAGE__)
endif ()


set( mri_image_header_files
        mri_image_export.h
        hoMRImage.h 
        hoMRImage.hxx )


add_library(gadgetron_toolbox_mri_image INTERFACE)


target_link_libraries(gadgetron_toolbox_mri_image
                    INTERFACE
                    gadgetron_toolbox_cpucore
                    gadgetron_toolbox_cpucore_math
                     )
target_include_directories(gadgetron_toolbox_mri_image
       INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

install(TARGETS gadgetron_toolbox_mri_image EXPORT gadgetron-export)


install(FILES ${mri_image_header_files}
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)