if (ARMADILLO_FOUND)
  install (FILES 
    spiral_flow_gpusense_cg.xml 
    spiral_flow_gpusense_sb.xml 
    spiral_interactive.xml 
    DESTINATION config)
elseif (ARMADILLO_FOUND)
  MESSAGE("Armadillo not found, only unoptimized spiral config files will be available")
endif (ARMADILLO_FOUND)

install (FILES 
  spiral_flow_gpusense_cg_unoptimized.xml 
  spiral_flow_gpusense_sb_unoptimized.xml 
  DESTINATION config)
