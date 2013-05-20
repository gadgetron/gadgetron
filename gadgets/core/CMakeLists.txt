IF (WIN32)
  ADD_DEFINITIONS(-DTIXML_USE_STL)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_CORE__)
  ADD_DEFINITIONS(-D__BUILD_GADGETS__)
ENDIF (WIN32)

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

find_package(Ismrmrd REQUIRED)
find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

add_library(gadgetron_core SHARED AcquisitionPassthroughGadget.cpp
  AcquisitionFinishGadget.cpp 
  AccumulatorGadget.cpp
  FFTGadget.cpp
  ImageFinishGadget.cpp
  CropAndCombineGadget.cpp
  ImageWriterGadget.cpp
  MRIImageWriter.cpp
  NoiseAdjustGadget.cpp
  ExtractGadget.cpp
  FloatToUShortGadget.cpp
  RemoveROOversamplingGadget.cpp
  CoilReductionGadget.cpp
  AutoScaleGadget.cpp
  FlowPhaseSubtractionGadget.cpp
  GadgetIsmrmrdReadWrite.cpp
  ${ISMRMRD_XSD_SOURCE}
  )

target_link_libraries(gadgetron_core 
  cpucore 
  ${Boost_LIBRARIES}
  ${ISMRMRD_LIBRARIES} ${FFTW3_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  ${XERCESC_LIBRARIES} 
  )

install (FILES	
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
  gadgetron_core_export.h
  DESTINATION include)

install(TARGETS gadgetron_core DESTINATION lib)
install(FILES default.xml default_short.xml DESTINATION config)

install(FILES ${ISMRMRD_SCHEMA_DIR}/ismrmrd.xsd DESTINATION schema)
install(FILES ${ISMRMRD_LIBRARIES} DESTINATION lib)
