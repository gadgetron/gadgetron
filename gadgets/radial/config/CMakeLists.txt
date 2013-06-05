if (ARMADILLO_FOUND)
  install (FILES 
    fixed_radial_mode0_realtime.xml 
    fixed_radial_mode1_realtime.xml 
    golden_radial_mode2_realtime.xml 
#    fixed_radial_gpusense_cg.xml 
#    fixed_radial_gpusense_sb.xml 
#    fixed_radial_gpusense_gp.xml
#    golden_radial_gpusense_cg.xml 
#    golden_radial_gpusense_sb.xml 
#    golden_radial_gpusense_gp.xml
    DESTINATION config)
elseif (ARMADILLO_FOUND)
  MESSAGE("Armadillo not found, no Blas/Lapack optimized config files are being installed!")
endif (ARMADILLO_FOUND)

install (FILES 
#  fixed_radial_gpusense_cg_unoptimized.xml 
#  fixed_radial_gpusense_sb_unoptimized.xml 
#  fixed_radial_gpusense_gp_unoptimized.xml
  DESTINATION config)
