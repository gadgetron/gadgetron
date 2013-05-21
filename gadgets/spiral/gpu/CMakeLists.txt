IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUSPIRAL__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)
find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

include_directories(
  ${ISMRMRD_XSD_INCLUDE_DIR}
  ${CUDA_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/gpu
  ${CMAKE_SOURCE_DIR}/gadgets/cgsense
  ${CMAKE_SOURCE_DIR}/gadgets/core
  ${CMAKE_SOURCE_DIR}/gadgets/sense
  ${CMAKE_SOURCE_DIR}/gadgets/spiral/cpu
  )

add_library(gadgetron_gpuspiral SHARED 
  gpuSpiralSensePrepGadget.cpp 
  ${ISMRMRD_XSD_SOURCE})

target_link_libraries(gadgetron_gpuspiral
  gpunfft gpusolvers gpuoperators gpuparallelmri gadgetron_cpuspiral
  ${ISMRMRD_LIBRARIES} ${XERCESC_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

install (TARGETS gadgetron_gpuspiral DESTINATION lib)

if (ARMADILLO_FOUND)
  install (FILES 
    spiral_gpusense_cg.xml 
    spiral_gpusense_sb.xml 
    spiral_gpusense_gp.xml
    DESTINATION config)
elseif (ARMADILLO_FOUND)
  MESSAGE("Armadillo not found, hence Gadgetron's Blas/Lapack optimized config files are not being installed")
endif (ARMADILLO_FOUND)

install (FILES 
  spiral_gpusense_cg_unoptimized.xml 
  spiral_gpusense_sb_unoptimized.xml 
  spiral_gpusense_gp_unoptimized.xml
  DESTINATION config)
