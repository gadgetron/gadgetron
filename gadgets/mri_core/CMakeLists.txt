IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_MRICORE__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)
find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${ARMADILLO_INCLUDE_DIRS}
)

add_library(gadgetron_mricore SHARED 
  MRINoiseAdjustGadget.cpp
  PCACoilGadget.cpp
  ${ISMRMRD_XSD_SOURCE}
  )

target_link_libraries(gadgetron_mricore cpucore cpucore_math 
  ${ISMRMRD_LIBRARIES} 
  ${FFTW3_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
  ${ARMADILLO_LIBRARIES}
  ${BOOST_LIBRARIES}
  ${XERCESC_LIBRARIES} 
  )

install(FILES MRINoiseAdjustGadget.h PCACoilGadget.h DESTINATION include)
install(FILES default_optimized.xml DESTINATION config)
install(TARGETS gadgetron_mricore DESTINATION lib)
