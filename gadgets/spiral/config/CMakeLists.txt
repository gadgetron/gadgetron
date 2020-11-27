if (CUDA_FOUND)
    if (ARMADILLO_FOUND)
      install (FILES 
        spiral_flow_gpusense_cg.xml 
        spiral_flow_gpusense_sb.xml 
        spiral_flow_generic_gpusense_cg.xml 
        spiral_flow_generic_gpusense_sb.xml 
        spiral_interactive.xml
	deblurring_recon_acctrig.xml
        DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
    elseif (ARMADILLO_FOUND)
      message("Armadillo not found, only unoptimized spiral config files will be available")
    endif ()

    install (FILES 
      spiral_flow_gpusense_cg_unoptimized.xml 
      spiral_flow_gpusense_sb_unoptimized.xml 
      DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
endif()
