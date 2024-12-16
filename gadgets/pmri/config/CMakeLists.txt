if (ARMADILLO_FOUND)
  install (FILES 
    generic_gpusense_cg.xml 
    generic_gpusense_cg_singleshot.xml 
    generic_gpusense_sb_singleshot.xml 
    generic_gpusense_nlcg_singleshot.xml
    generic_gpu_ktsense_singleshot.xml 
    DESTINATION ${PINGVIN_INSTALL_CONFIG_PATH} COMPONENT main)
elseif (ARMADILLO_FOUND)
  message("Armadillo not found, only unoptimized generic trajectory config files will be available")
endif ()

#install (FILES 
#  generic_gpusense_cg_unoptimized.xml 
#  generic_gpusense_sb_unoptimized.xml 
#  DESTINATION ${PINGVIN_INSTALL_CONFIG_PATH})
