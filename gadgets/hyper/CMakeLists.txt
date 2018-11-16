if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_HYPER__)
endif ()

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/hyper
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/dwt/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/gadgets/pmri
  ${HDF5_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIR}/cpp
  ${ARMADILLO_INCLUDE_DIRS}
)

add_library(gadgetron_hyper SHARED 
  NFFT2DGadget.h NFFT2DGadget.cpp
  CMRT3DGadget.h CMRT3DGadget.cpp
  CMRTGadget.h CMRTGadget.cpp
  CSIGadget.h CSIGadget.cpp
  gpuCSICoilEstimationGadget.cpp
  
  )

set_target_properties(gadgetron_hyper PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_hyper
  gadgetron_gadgetbase
  gadgetron_toolbox_gpunfft
  gadgetron_toolbox_gpucore
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_cpufft
  gadgetron_toolbox_gpudwt
  gadgetron_toolbox_cpucore_math
  gadgetron_toolbox_gpusolvers
  gadgetron_toolbox_gpuoperators
  gadgetron_toolbox_gpuparallelmri
  gadgetron_toolbox_hyper
  ${ISMRMRD_LIBRARIES}  
  ${FFTW3_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
  ${BOOST_LIBRARIES}
  )

install(FILES 
  NFFT2DGadget.h
  CMRTGadget.h
  CMRT3DGadget.h
  CSIGadget.h
  gadgetron_hyper_export.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH})

install(FILES NFFT2D.xml CMRT3D.xml CMRT.xml FS-CSI.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH})

install(TARGETS gadgetron_hyper DESTINATION lib COMPONENT main)
