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
    fixed_radial_mode0_gpu_ktsense.xml 
    fixed_radial_mode1_gpu_ktsense.xml 
    golden_radial_mode2_gpu_ktsense.xml 
    DESTINATION config)
elseif (ARMADILLO_FOUND)
  MESSAGE("Armadillo not found, only unoptimized radial config files will be available")
endif (ARMADILLO_FOUND)

install (FILES 
  fixed_radial_mode0_gpusense_cg_unoptimized.xml 
  fixed_radial_mode1_gpusense_cg_unoptimized.xml 
  golden_radial_mode2_gpusense_cg_unoptimized.xml 
  fixed_radial_mode0_gpusense_sb_unoptimized.xml 
  fixed_radial_mode1_gpusense_sb_unoptimized.xml 
  golden_radial_mode2_gpusense_sb_unoptimized.xml 
  DESTINATION config)
