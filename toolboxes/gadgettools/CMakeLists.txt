include_directories(${ACE_INCLUDE_DIR}  
			       ${CMAKE_SOURCE_DIR}/apps/gadgetron)

add_library(gadgettools SHARED GadgetronConnector.cpp)

target_link_libraries(gadgettools ${ACE_LIBRARIES})

install(TARGETS gadgettools DESTINATION lib)

install (FILES	GadgetSocketReceiver.h 
			GadgetSocketSender.h
			GadgetXml.h
			GadgetronConnector.h
			DESTINATION include)
