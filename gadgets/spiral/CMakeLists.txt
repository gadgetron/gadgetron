IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUSPIRAL__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)
find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/gpu
  ${CMAKE_SOURCE_DIR}/gadgets/cgsense
  ${CMAKE_SOURCE_DIR}/gadgets/core
  ${CMAKE_SOURCE_DIR}/gadgets/sense
  ${ISMRMRD_XSD_INCLUDE_DIR}
  ${CUDA_INCLUDE_DIRS}
  )

add_library(gadgetron_gpuspiral SHARED 
  vds.cpp 
  gpuSpiralSensePrepGadget.cpp 
  ${ISMRMRD_XSD_SOURCE})

target_link_libraries(gadgetron_gpuspiral
  cpucore gpucore gpunfft gpusolvers gpuoperators gpuparallelmri
  ${ISMRMRD_LIBRARIES} ${XERCESC_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

install (TARGETS gadgetron_gpuspiral DESTINATION lib)
install (FILES vds.h DESTINATION include)

if (ARMADILLO_FOUND)
  install (FILES 
    spiral_gpusense_cg.xml 
    spiral_gpusense_sb.xml 
    DESTINATION config)
elseif (ARMADILLO_FOUND)
  MESSAGE("Armadillo not found. Hence only a limited number (of unoptimized) config files will be available!")
endif (ARMADILLO_FOUND)

install (FILES 
  spiral_gpusense_cg_unoptimized.xml 
  spiral_gpusense_sb_unoptimized.xml 
  DESTINATION config)
