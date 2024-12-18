add_library(pingvin_epi SHARED
  EPIReconXGadget.h EPIReconXGadget.cpp
  EPICorrGadget.h EPICorrGadget.cpp
  EPIPackNavigatorGadget.h EPIPackNavigatorGadget.cpp
  FFTXGadget.h FFTXGadget.cpp
  CutXGadget.h CutXGadget.cpp
  OneEncodingGadget.h OneEncodingGadget.cpp
  epi.xml
  epi_gtplus_grappa.xml
)

set_target_properties(pingvin_epi PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(
  pingvin_epi
  pingvin_core
  pingvin_mricore
  pingvin_toolbox_cpucore
  pingvin_toolbox_cpufft
  pingvin_toolbox_cpucore_math
  pingvin_toolbox_log
  pingvin_toolbox_epi
)


install(FILES
  EPIReconXGadget.h
  EPICorrGadget.h
  EPIPackNavigatorGadget.h
  FFTXGadget.h
  DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS pingvin_epi
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
  epi.xml
  epi_gtplus_grappa.xml
  DESTINATION ${PINGVIN_INSTALL_CONFIG_PATH} COMPONENT main)