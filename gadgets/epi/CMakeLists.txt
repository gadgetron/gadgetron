if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_EPI__)
endif ()


add_library(gadgetron_epi SHARED 
  EPIReconXGadget.h EPIReconXGadget.cpp
  EPICorrGadget.h EPICorrGadget.cpp
  EPIPackNavigatorGadget.h EPIPackNavigatorGadget.cpp
  FFTXGadget.h FFTXGadget.cpp
  CutXGadget.h CutXGadget.cpp
  OneEncodingGadget.h OneEncodingGadget.cpp
  epi.xml
  epi_gtplus_grappa.xml
)

set_target_properties(gadgetron_epi PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(
  gadgetron_epi
  gadgetron_core
  gadgetron_mricore
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_cpufft
  gadgetron_toolbox_cpucore_math
  gadgetron_toolbox_log
  gadgetron_toolbox_epi
  ${ISMRMRD_LIBRARIES}
)


install(FILES
  EPIReconXGadget.h
  EPICorrGadget.h
  EPIPackNavigatorGadget.h
  FFTXGadget.h
  gadgetron_epi_export.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_epi
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
  epi.xml
  epi_gtplus_grappa.xml
  DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)