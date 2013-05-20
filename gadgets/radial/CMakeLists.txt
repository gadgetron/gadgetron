#add_subdirectory(cpu)
  
if (CUDA_FOUND)
  add_subdirectory(gpu)
endif(CUDA_FOUND)

install(FILES
  radial_single_optimized.xml 
  radial_single_fixed_optimized.xml
  radial_ms_optimized.xml
  radial_sb.xml
  radial_gp.xml
  DESTINATION config)
