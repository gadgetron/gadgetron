include_directories(${ACE_INCLUDE_DIR} ${TINYXML_INCLUDE_DIR})

add_executable(gadgetron main.cpp GadgetStreamController.cpp GadgetServerAcceptor.cpp)

target_link_libraries(gadgetron ${ACE_LIBRARIES} ${TINYXML_LIBRARIES})

install(TARGETS gadgetron DESTINATION bin)

install(FILES 	
	GadgetContainerMessage.h
	Gadget.h
	GadgetServerAcceptor.h
	GadgetStreamController.h
	Gadgetron.h
	GadgetMessageInterface.h
DESTINATION include)

