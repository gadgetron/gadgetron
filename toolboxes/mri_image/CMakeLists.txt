if (WIN32)
    ADD_DEFINITIONS(-D__BUILD_GADGETRON_MRI_IMAGE__)
endif (WIN32)

include_directories(
    ${Boost_INCLUDE_DIR} 
    ${ISMRMRD_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_image
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/klt/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
)

set( mri_image_header_files
        mri_image_export.h
        hoMRImage.h 
        hoMRImage.hxx )

set( mri_image_source_files )

#add_library(gadgetron_toolbox_mri_image SHARED ${mri_image_header_files} ${mri_image_source_files} )

#set_target_properties(gadgetron_toolbox_mri_image PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
#set_target_properties(gadgetron_toolbox_mri_image PROPERTIES LINKER_LANGUAGE CXX)

#target_link_libraries(gadgetron_toolbox_mri_image 
#                    gadgetron_toolbox_cpucore 
#                    gadgetron_toolbox_cpucore_math 
#                    ${ARMADILLO_LIBRARIES} )

# install(TARGETS gadgetron_toolbox_mri_image DESTINATION lib COMPONENT main)

install(FILES ${mri_image_header_files}
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
