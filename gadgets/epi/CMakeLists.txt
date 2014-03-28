IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_EPI__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)
find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/epi
  ${ISMRMRD_XSD_INCLUDE_DIR}
  ${ARMADILLO_INCLUDE_DIRS}
)

add_library(gadgetron_epi SHARED 
  EPIReconXGadget.h EPIReconXGadget.cpp
  ${ISMRMRD_XSD_SOURCE}
)

target_link_libraries(
  gadgetron_epi
  cpucore
  cpucore_math
  ${ISMRMRD_LIBRARIES} 
  ${XERCESC_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
)

install(FILES 
  EPIReconXGadget.h
  gadgetron_epi_export.h
  DESTINATION include)

install(TARGETS gadgetron_epi DESTINATION lib)

install(FILES
  epi.xml
  DESTINATION config)

