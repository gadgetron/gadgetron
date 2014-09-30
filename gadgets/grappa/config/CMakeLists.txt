install (FILES grappa_unoptimized.xml grappa_unoptimized_float.xml DESTINATION config COMPONENT main)

if(ARMADILLO_FOUND)
  install (FILES grappa.xml grappa_float.xml grappa_device.xml DESTINATION config COMPONENT main)
endif(ARMADILLO_FOUND)
