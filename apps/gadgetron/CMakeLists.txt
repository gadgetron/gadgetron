IF (WIN32)
ADD_DEFINITIONS(-DTIXML_USE_STL)
ENDIF (WIN32)

include_directories(
  ${ACE_INCLUDE_DIR} 
  ${Boost_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/apps/gadgetron
  ${CMAKE_SOURCE_DIR}/dependencies/tinyxml
  ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools)

add_executable(gadgetron main.cpp GadgetStreamController.cpp GadgetServerAcceptor.cpp)

target_link_libraries(gadgetron tinyxml optimized gadgettools ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

install(TARGETS gadgetron DESTINATION bin)

install(FILES gadgetron_start.pl DESTINATION bin PERMISSIONS OWNER_WRITE OWNER_EXECUTE OWNER_READ GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install(FILES 	
	GadgetContainerMessage.h
	Gadget.h
	GadgetServerAcceptor.h
	GadgetStreamController.h
	Gadgetron.h
	GadgetMessageInterface.h
	EndGadget.h
	GadgetronExport.h
DESTINATION include)
