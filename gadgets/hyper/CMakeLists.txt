IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_HYPER__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/hyper
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${HDF5_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIR}/cpp
  ${ARMADILLO_INCLUDE_DIRS}
)

add_library(gadgetron_hyper SHARED 
  NFFT2DGadget.h NFFT2DGadget.cpp
  CMRT3DGadget.h CMRT3DGadget.cpp
  CMRTGadget.h CMRTGadget.cpp
  )

target_link_libraries(gadgetron_hyper 
  gpunfft
  gpucore
  cpucore
  cpucore_math
  gpusolvers
  ${ISMRMRD_LIBRARIES} 
  ${FFTW3_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
  ${BOOST_LIBRARIES}
  )

install(FILES 
  NFFT2DGadget.h 
  gadgetron_hyper_export.h
  DESTINATION include)

install(FILES NFFT2D.xml CMRT3D.xml CMRT.xml DESTINATION config)

install(TARGETS gadgetron_hyper DESTINATION lib)
