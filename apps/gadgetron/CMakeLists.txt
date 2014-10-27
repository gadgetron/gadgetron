configure_file(gadgetron_config.in gadgetron_config.h)

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/apps/gadgetron
  ${CMAKE_SOURCE_DIR}/toolboxes/cloudbus
  ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
  ${Boost_INCLUDE_DIR}
  ${ACE_INCLUDE_DIR}
  )

add_executable(gadgetron 
  main.cpp 
  GadgetServerAcceptor.h
  GadgetServerAcceptor.cpp 
  GadgetStreamController.h
  EndGadget.h 
  Gadget.h 
  GadgetContainerMessage.h 
  GadgetMessageInterface.h 
  Gadgetron.h 
  GadgetronExport.h 
  gadgetron_xml.h
  gadgetron_xml.cpp
  pugixml.cpp
  )

target_link_libraries(gadgetron 
  gadgetron_gadgetbase
  gadgetron_toolbox_gadgettools gadgetron_toolbox_cloudbus 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
 )

add_executable(gadgetron_info
  gadgetron_info.cpp
)

target_link_libraries(gadgetron_info 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
 )

add_library(gadgetron_gadgetbase SHARED
  Gadget.cpp
  GadgetStreamController.cpp
)

target_link_libraries(gadgetron_gadgetbase
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  gadgetron_toolbox_gadgettools
)

set_target_properties (gadgetron_gadgetbase PROPERTIES COMPILE_DEFINITIONS "__BUILD_GADGETRON_GADGETBASE__")

install(TARGETS gadgetron gadgetron_info DESTINATION bin COMPONENT main)
install(TARGETS gadgetron_gadgetbase DESTINATION lib COMPONENT main)

install(FILES 
  EndGadget.h
  Gadget.h
  GadgetContainerMessage.h
  GadgetMessageInterface.h
  Gadgetron.h
  GadgetronExport.h
  gadgetron_paths.h
  gadgetron_xml.h
  GadgetServerAcceptor.h
  GadgetStreamController.h
  ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_config.h
  DESTINATION include COMPONENT main) 

install(FILES 
  gadgetron.xml.example
  DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

install(FILES 
  schema/gadgetron.xsd
  DESTINATION schema COMPONENT main)

add_subdirectory(webapp)
