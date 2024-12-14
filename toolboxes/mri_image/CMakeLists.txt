set( mri_image_header_files
        hoMRImage.h
        hoMRImage.hxx )

add_library(pingvin_toolbox_mri_image INTERFACE)

target_link_libraries(pingvin_toolbox_mri_image
                    INTERFACE
                    pingvin_toolbox_cpucore
                    pingvin_toolbox_cpucore_math
                     )
target_include_directories(pingvin_toolbox_mri_image
       INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

install(TARGETS pingvin_toolbox_mri_image EXPORT pingvin-export)


install(FILES ${mri_image_header_files}
  DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)