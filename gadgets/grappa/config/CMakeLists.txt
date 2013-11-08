install (FILES grappa_unoptimized.xml grappa_unoptimized_float.xml DESTINATION config)

if(ARMADILLO_FOUND)
  install (FILES grappa.xml grappa_float.xml DESTINATION config)
endif(ARMADILLO_FOUND)
