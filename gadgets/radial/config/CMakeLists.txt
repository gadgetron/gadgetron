if (ARMADILLO_FOUND)
  install (FILES 
    fixed_radial_mode0_realtime.xml 
    fixed_radial_mode1_realtime.xml 
    golden_radial_mode2_realtime.xml 
    fixed_radial_mode0_gpusense_cg.xml 
    fixed_radial_mode1_gpusense_cg.xml 
    golden_radial_mode2_gpusense_cg.xml 
    fixed_radial_mode0_gpusense_sb.xml 
    fixed_radial_mode1_gpusense_sb.xml 
    golden_radial_mode2_gpusense_sb.xml 
    DESTINATION config)
elseif (ARMADILLO_FOUND)
  MESSAGE("Armadillo not found. Hence only a limited number (of unoptimized) config files will be available!")
endif (ARMADILLO_FOUND)

install (FILES 
  fixed_radial_mode0_gpusense_cg_unoptimized.xml 
  fixed_radial_mode1_gpusense_cg_unoptimized.xml 
  golden_radial_mode2_gpusense_cg_unoptimized.xml 
  fixed_radial_mode0_gpusense_sb_unoptimized.xml 
  fixed_radial_mode1_gpusense_sb_unoptimized.xml 
  golden_radial_mode2_gpusense_sb_unoptimized.xml 
  DESTINATION config)
