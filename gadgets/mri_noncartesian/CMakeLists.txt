if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_MRI_NONCARTESIAN__)
endif ()

set( gadgetron_mri_noncartesian_header_files
	CPUGriddingReconGadget.h NonCartesianTools.h GriddingReconGadgetBase.h GriddingReconGadgetBase.hpp)

set( gadgetron_mri_noncartesian_src_files
	CPUGriddingReconGadget.cpp NonCartesianTools.cpp)

set( gadgetron_mri_noncartesian_config_files
	config/Generic_CPU_Gridding_Recon.xml
		config/Generic_Spiral.xml
		config/Generic_Spiral_SNR.xml
		config/Generic_Spiral_Flag.xml
 )


if(CUDA_FOUND)
	list(APPEND gadgetron_mri_noncartesian_header_files GriddingReconGadget.h)
	list(APPEND gadgetron_mri_noncartesian_src_files GriddingReconGadget.cpp)

	list(APPEND gadgetron_mri_noncartesian_config_files config/Generic_Spiral.xml config/Generic_Spiral_SNR.xml)

endif()

add_library(gadgetron_mri_noncartesian SHARED 
  gadgetron_mri_noncartesian_export.h 
  ${gadgetron_mri_noncartesian_header_files} 
  ${gadgetron_mri_noncartesian_src_files}
  ${gadgetron_mri_noncartesian_config_files}
)

set_target_properties(gadgetron_mri_noncartesian PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_mri_noncartesian 
    gadgetron_core
    gadgetron_toolbox_log
    gadgetron_toolbox_mri_core
    gadgetron_mricore
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_cpufft
    gadgetron_toolbox_cpunfft
    gadgetron_toolbox_mri_core
    gadgetron_toolbox_cpuoperator
    gadgetron_toolbox_image_analyze_io


)
if(CUDA_FOUND)
    target_link_libraries(gadgetron_mri_noncartesian
	gadgetron_toolbox_gpu
	)
endif()

install(FILES 
    gadgetron_mri_noncartesian_export.h
    ${gadgetron_mri_noncartesian_header_files}
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_mri_noncartesian
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES ${gadgetron_mri_noncartesian_config_files} DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)