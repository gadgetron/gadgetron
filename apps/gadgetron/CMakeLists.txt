include_directories(
  ${CMAKE_SOURCE_DIR}/apps/gadgetron
  ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
  ${Boost_INCLUDE_DIR}
  ${ACE_INCLUDE_DIR} 
  ${XSD_INCLUDE_DIR}
  )

#Process the XSD files
SET(XSDS schema/gadgetron.xsd)
SET(XSD_ARGS cxx-tree --generate-serialization)
WRAP_XSD(XSDS_SOURCES XSD_INCLUDES ${CMAKE_CURRENT_BINARY_DIR}/schema ${XSDS} OPTIONS ${XSD_ARGS})
INCLUDE_DIRECTORIES(${XSD_INCLUDES} ${XERCESC_INCLUDE_DIR})

add_executable(gadgetron 
  ${XSDS_SOURCES} 
  main.cpp 
  EndGadget.h
  Gadget.h
  GadgetContainerMessage.h
  GadgetMessageInterface.h
  Gadgetron.h
  GadgetronExport.h
  )

target_link_libraries(gadgetron 
  gadgettools 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
  ${XERCESC_LIBRARIES}
 )

install(TARGETS gadgetron DESTINATION bin)

install(FILES 	
  GadgetContainerMessage.h
  Gadget.h
  Gadgetron.h
  GadgetMessageInterface.h
  EndGadget.h
  GadgetronExport.h
  DESTINATION include) 

install(FILES 	
  gadgetron.xml.example
  DESTINATION config)

install(FILES 	
  schema/gadgetron.xsd
  DESTINATION schema)
