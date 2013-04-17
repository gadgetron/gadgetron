IF (WIN32)
  ADD_DEFINITIONS(-DTIXML_USE_STL)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_MRICORE__)
  ADD_DEFINITIONS(-D__BUILD_GADGETS__)
ENDIF (WIN32)

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

find_package(Ismrmrd REQUIRED)
find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${ARMADILLO_INCLUDE_DIRS}
)

add_library(gadgetronmricore SHARED 
  MRINoiseAdjustGadget.cpp
  PCACoilGadget.cpp
  ${ISMRMRD_XSD_SOURCE}
  )

target_link_libraries(gadgetronmricore cpucore cpucore_math 
  ${ISMRMRD_LIBRARIES}  
  ${FFTW3_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
  ${ARMADILLO_LIBRARIES}
  ${BOOST_LIBRARIES}
  ${XERCESC_LIBRARIES} 
  )

install(FILES MRINoiseAdjustGadget.h PCACoilGadget.h DESTINATION include)
install(FILES default_optimized.xml radial_single_optimized.xml radial_single_fixed_optimized.xml radial_ms_optimized.xml DESTINATION config)
install(TARGETS gadgetronmricore DESTINATION lib)
