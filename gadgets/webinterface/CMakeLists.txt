IF (WIN32)
ADD_DEFINITIONS(-DTIXML_USE_STL)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_WEBINTERFACE__)
ADD_DEFINITIONS(-D__BUILD_GADGETS__)
ENDIF (WIN32)

include_directories(   ${ACE_INCLUDE_DIR} 
				       ${Boost_INCLUDE_DIR}
					   ${CMAKE_SOURCE_DIR}/dependencies/tinyxml
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron
					   ${CMAKE_SOURCE_DIR}/toolboxes/hostutils
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools)

add_library(gadgetronwebinterface SHARED WebInterfaceGadget.cpp
								 mongoose.c)

target_link_libraries(gadgetronwebinterface tinyxml optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

install (FILES	WebInterfaceGadget.h
			gadgetronwebinterface_export.h
			DESTINATION include)
			
install (FILES	gca.html
			DESTINATION html)
			
install (FILES	radial_ms_webinterface.xml
			DESTINATION config)
			

install(TARGETS gadgetronwebinterface DESTINATION lib)
