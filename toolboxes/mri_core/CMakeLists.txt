if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_MRI_CORE__)
endif ()



set( mri_core_header_files
        mri_core_export.h
        mri_core_def.h
        mri_core_data.h
        mri_core_utility.h
        mri_core_kspace_filter.h
        mri_core_grappa.h 
        mri_core_spirit.h 
        mri_core_coil_map_estimation.h 
        mri_core_dependencies.h 
        mri_core_acquisition_bucket.h
        mri_core_girf_correction.h
        mri_core_partial_fourier.h )

set( mri_core_source_files
        mri_core_utility.cpp 
        mri_core_grappa.cpp 
        mri_core_spirit.cpp 
        mri_core_kspace_filter.cpp
        mri_core_coil_map_estimation.cpp 
        mri_core_dependencies.cpp
        mri_core_girf_correction.cpp
        mri_core_partial_fourier.cpp )

add_library(gadgetron_toolbox_mri_core SHARED 
     ${mri_core_header_files} ${mri_core_source_files} )

set_target_properties(gadgetron_toolbox_mri_core PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
set_target_properties(gadgetron_toolbox_mri_core PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(gadgetron_toolbox_mri_core
                    gadgetron_core
                    gadgetron_toolbox_cpucore 
                    gadgetron_toolbox_cpucore_math 
                    armadillo
                    boost
                    gadgetron_toolbox_cpufft 
                    gadgetron_toolbox_cpuklt )
target_include_directories(
        gadgetron_toolbox_mri_core
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${ISMRMRD_INCLUDE_DIR}
)
install(TARGETS gadgetron_toolbox_mri_core DESTINATION lib COMPONENT main)

install(FILES
  ${mri_core_header_files}
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
