if(CUDA_FOUND)
  if(CPU_REG)
    install (FILES 
      cpureg_cartesian_averaging.xml
      DESTINATION config COMPONENT main)
  endif(CPU_REG)
  
  if(GPU_REG)
    install (FILES 
      gpureg_cartesian_averaging.xml
      DESTINATION config COMPONENT main)
  endif(GPU_REG)
endif(CUDA_FOUND)
