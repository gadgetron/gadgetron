include_directories(${ACE_INCLUDE_DIR} 
				    ${Boost_INCLUDE_DIR}  
			        ${CMAKE_SOURCE_DIR}/apps/gadgetron)

add_library(gadgettools SHARED GadgetronConnector.cpp)

target_link_libraries(gadgettools optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

install(TARGETS gadgettools DESTINATION lib)

install (FILES	GadgetSocketReceiver.h 
			GadgetSocketSender.h
			GadgetXml.h
			GadgetronConnector.h
			gadgettools_export.h
			GadgetronSlotContainer.h
			DESTINATION include)
