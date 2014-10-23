
install (FILES grappa_unoptimized.xml grappa_unoptimized_float.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

if(ARMADILLO_FOUND)
  install (FILES grappa.xml grappa_float.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
endif(ARMADILLO_FOUND)
