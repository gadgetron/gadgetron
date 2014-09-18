IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_EPI__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/epi
  ${ARMADILLO_INCLUDE_DIRS}
)

add_library(gadgetron_epi SHARED 
  EPIReconXGadget.h EPIReconXGadget.cpp
  EPICorrGadget.h EPICorrGadget.cpp
  FFTXGadget.h FFTXGadget.cpp
  CutXGadget.h CutXGadget.cpp
  epi.xml
  epi_gtplus_grappa.xml
)

target_link_libraries(
  gadgetron_epi
  cpucore
  cpufft
  cpucore_math
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
)

install(FILES 
  EPIReconXGadget.h
  EPICorrGadget.h
  FFTXGadget.h
  gadgetron_epi_export.h
  DESTINATION include)

install(TARGETS gadgetron_epi DESTINATION lib)

install(FILES
  epi.xml
  epi_gtplus_grappa.xml
  DESTINATION config)

