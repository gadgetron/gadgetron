include_directories(
  ${CMAKE_SOURCE_DIR}/apps/gadgetron
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
  GadgetStreamController.cpp 
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
  gadgettools 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
 )

install(TARGETS gadgetron DESTINATION bin)

install(FILES 
  GadgetContainerMessage.h
  Gadget.h
  Gadgetron.h
  GadgetMessageInterface.h
  EndGadget.h
  GadgetServerAcceptor.h
  GadgetStreamController.h
  GadgetronExport.h
  DESTINATION include) 

install(FILES 
  gadgetron.xml.example
  DESTINATION config)

install(FILES 
  schema/gadgetron.xsd
  DESTINATION schema)
