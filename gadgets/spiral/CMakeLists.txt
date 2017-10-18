if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_SPIRAL__)
endif ()

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/gadgets/pmri
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/gpu
  ${CUDA_INCLUDE_DIRS}
  )

set(SPIRAL_HEADER_FILES gadgetron_spiral_export.h SpiralToGenericGadget.h vds.h)
set(SPIRAL_SRC_FILES vds.cpp  SpiralToGenericGadget.cpp)

if (CUDA_FOUND)
    set(SPIRAL_HEADER_FILES ${SPIRAL_HEADER_FILES} gpuSpiralSensePrepGadget.h)
    set(SPIRAL_SRC_FILES ${SPIRAL_SRC_FILES}  gpuSpiralSensePrepGadget.cpp)
endif()
  
add_library(gadgetron_spiral SHARED ${SPIRAL_HEADER_FILES} ${SPIRAL_SRC_FILES})

set_target_properties(gadgetron_spiral PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_spiral
  gadgetron_gadgetbase
  gadgetron_toolbox_log
  gadgetron_toolbox_cpucore 
  ${ISMRMRD_LIBRARIES} 
  ${FFTW3_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

if (CUDA_FOUND)
    target_link_libraries(gadgetron_spiral 
                                ${CUDA_LIBRARIES} 
                                gadgetron_toolbox_gpucore 
                                gadgetron_toolbox_gpunfft 
                                gadgetron_toolbox_gpusolvers 
                                gadgetron_toolbox_gpuoperators
                                gadgetron_toolbox_gpuparallelmri )
endif()

install (TARGETS gadgetron_spiral DESTINATION lib COMPONENT main)
install (FILES ${SPIRAL_HEADER_FILES}
                     DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

add_subdirectory(config)
