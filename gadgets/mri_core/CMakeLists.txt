IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_MRICORE__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)
find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${ARMADILLO_INCLUDE_DIRS}
)

if (ARMADILLO_FOUND)
  list(APPEND OPTIMIZED_GADGETS NoiseAdjustGadget.cpp)
  list(APPEND OPTIMIZED_GADGETS PCACoilGadget.cpp)
  list(APPEND OPTIMIZED_GADGET_HEADERS NoiseAdjustGadget.h)
  list(APPEND OPTIMIZED_GADGET_HEADERS PCACoilGadget.h)
  list(APPEND OPTIMIZED_GADGET_LIBS cpucore_math ${ARMADILLO_LIBRARIES})
elseif (ARMADILLO_FOUND)
  MESSAGE("Armadillo not found, i.e. not compiling Blas/Lapack optimized MRI Gadgets")
endif (ARMADILLO_FOUND)

add_library(gadgetron_mricore SHARED 
  AcquisitionPassthroughGadget.cpp
  AcquisitionFinishGadget.cpp 
  AccumulatorGadget.cpp
  FFTGadget.cpp
  ImageFinishGadget.cpp
  CropAndCombineGadget.cpp
  ImageWriterGadget.cpp
  MRIImageWriter.cpp
  NoiseAdjustGadget_unoptimized.cpp
  ExtractGadget.cpp
  FloatToUShortGadget.cpp
  RemoveROOversamplingGadget.cpp
  CoilReductionGadget.cpp
  AutoScaleGadget.cpp
  FlowPhaseSubtractionGadget.cpp
  GadgetIsmrmrdReadWrite.cpp
  PhysioInterpolationGadget.cpp
  ${OPTIMIZED_GADGETS}
  ${ISMRMRD_XSD_SOURCE}
  )

target_link_libraries(gadgetron_mricore 
  cpucore
  ${ISMRMRD_LIBRARIES} 
  ${FFTW3_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
  ${BOOST_LIBRARIES}
  ${XERCESC_LIBRARIES} 
  ${OPTIMIZED_GADGET_LIBS}
  )

install(FILES 
  AccumulatorGadget.h
  AcquisitionFinishGadget.h
  AcquisitionPassthroughGadget.h
  CropAndCombineGadget.h	
  ExtractGadget.h
  FloatToUShortGadget.h
  FFTGadget.h
  GadgetMRIHeaders.h
  ImageFinishGadget.h
  ImageWriterGadget.h
  MRIImageWriter.h
  NoiseAdjustGadget.h
  RemoveROOversamplingGadget.h
  CoilReductionGadget.h
  AutoScaleGadget.h
  FlowPhaseSubtractionGadget.h
  GadgetIsmrmrdReadWrite.h
  gadgetron_mricore_export.h
  ${OPTIMIZED_GADGET_HEADERS}
  DESTINATION include)

if (ARMADILLO_FOUND)
  install(FILES default.xml default_short.xml default_optimized.xml DESTINATION config)
endif (ARMADILLO_FOUND)

install(TARGETS gadgetron_mricore DESTINATION lib)
install(FILES ${ISMRMRD_SCHEMA_DIR}/ismrmrd.xsd DESTINATION schema)
