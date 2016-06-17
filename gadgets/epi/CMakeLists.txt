IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_EPI__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
  ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/epi
  ${ARMADILLO_INCLUDE_DIRS}
)

add_library(gadgetron_epi SHARED 
  EPIReconXGadget.h EPIReconXGadget.cpp
  EPICorrGadget.h EPICorrGadget.cpp
  FFTXGadget.h FFTXGadget.cpp
  CutXGadget.h CutXGadget.cpp
  OneEncodingGadget.h OneEncodingGadget.cpp
  epi.xml
  epi_gtplus_grappa.xml
)

set_target_properties(gadgetron_epi PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(
  gadgetron_epi
  gadgetron_gadgetbase
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_cpufft
  gadgetron_toolbox_cpucore_math
  gadgetron_toolbox_log
  ${ISMRMRD_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
)

install(FILES 
  EPIReconXGadget.h
  EPICorrGadget.h
  FFTXGadget.h
  gadgetron_epi_export.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_epi DESTINATION lib COMPONENT main)

install(FILES
  epi.xml
  epi_gtplus_grappa.xml
  DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
