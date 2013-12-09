IF (WIN32)
ADD_DEFINITIONS(-DTIXML_USE_STL)
ENDIF (WIN32)

find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

include_directories(
  ${Boost_INCLUDE_DIR}
  ${ACE_INCLUDE_DIR} 
  ${XSD_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/apps/gadgetron
  ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils)

#Process the XSD files
SET(XSDS schema/gadgetron.xsd)
SET(XSD_ARGS cxx-tree --generate-serialization)
WRAP_XSD(XSDS_SOURCES XSD_INCLUDES ${CMAKE_CURRENT_BINARY_DIR}/schema ${XSDS} OPTIONS ${XSD_ARGS})
INCLUDE_DIRECTORIES(${XSD_INCLUDES} ${XERCESC_INCLUDE_DIR})

add_executable(gadgetron 
  ${XSDS_SOURCES} 
  main.cpp 
  #GadgetStreamController.cpp 
  #GadgetServerAcceptor.cpp 
  EndGadget.h
  Gadget.h
  GadgetContainerMessage.h
  GadgetMessageInterface.h
  Gadgetron.h
  GadgetronExport.h
  #GadgetServerAcceptor.h
  #GadgetStreamController.h 
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
  #GadgetServerAcceptor.h
  #GadgetStreamController.h
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
