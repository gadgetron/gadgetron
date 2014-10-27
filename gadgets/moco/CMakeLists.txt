if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_MOCO__)
endif (WIN32)

find_package(Ismrmrd REQUIRED)

set( CPU_REG 0)
set( GPU_REG 0)

if(ARMADILLO_FOUND)
  if(ARMADILLO_VERSION_STRING VERSION_GREATER "3.819" )
    
    set(CPU_REG 1)
    list(APPEND CPU_GADGETS cpuRegistrationAveragingGadget.cpp)
    list(APPEND CPU_LIBS gadgetron_toolbox_cpucore_math gadgetron_toolbox_cpureg)
    
  elseif(ARMADILLO_VERSION_STRING VERSION_GREATER "3.819" )
    MESSAGE("Armadillo of at least version 3.820 not found, not compiling cpu-based registration gadgets")
  endif(ARMADILLO_VERSION_STRING VERSION_GREATER "3.819" )
elseif (ARMADILLO_FOUND)
  MESSAGE("Armadillo not found, not compiling cpu-based registration gadgets")
endif (ARMADILLO_FOUND)

if(CUDA_FOUND)
  set(GPU_REG 1)
  list(APPEND GPU_GADGETS gpuRegistrationAveragingGadget.cpp gpuRegistrationScatteringGadget.cpp)
  list(APPEND GPU_LIBS gadgetron_toolbox_gpucore gadgetron_toolbox_gpureg ${CUDA_LIBRARIES})
elseif (CUDA_FOUND)
  MESSAGE("Cuda not found, not compiling gpu-based registration gadgets")
endif (CUDA_FOUND)

include_directories(   
  ${ACE_INCLUDE_DIR} 
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/gadgets/moco
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/
  ${HDF5_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIR}/cpp
)

if(CPU_REG)
  include_directories(   
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
    ${ARMADILLO_INCLUDE_DIRS}
    )
endif(CPU_REG)

if(GPU_REG)
  include_directories(   
    ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/gpu
    ${CUDA_INCLUDE_DIRS}
    )
endif(GPU_REG)

if(CUDA_FOUND)
    add_library(gadgetron_moco SHARED
        cpuRegistrationAveragingGadget.h
        gadgetron_moco_export.h
        gpuRegistrationAveragingGadget.h
        gpuRegistrationScatteringGadget.h
        RegistrationAveragingGadget.h
        RegistrationScatteringGadget.h
        ${CPU_GADGETS}
        ${GPU_GADGETS}
      )

    set_target_properties(gadgetron_moco PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})                                                                                                                                                                                                         

    target_link_libraries(gadgetron_moco
      gadgetron_gadgetbase
      gadgetron_toolbox_cpucore gadgetron_mricore ${CPU_LIBS} ${GPU_LIBS}
      ${Boost_LIBRARIES} ${ISMRMRD_LIBRARIES}
      optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
      )

    install (TARGETS gadgetron_moco DESTINATION lib COMPONENT main)
endif(CUDA_FOUND)

install(FILES   cpuRegistrationAveragingGadget.h
                gadgetron_moco_export.h
                gpuRegistrationAveragingGadget.h
                gpuRegistrationScatteringGadget.h
                RegistrationAveragingGadget.h
                RegistrationScatteringGadget.h
                DESTINATION include COMPONENT main)

add_subdirectory(config)
