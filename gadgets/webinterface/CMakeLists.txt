IF (WIN32)
ADD_DEFINITIONS(-DTIXML_USE_STL)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_WEBINTERFACE__)
ADD_DEFINITIONS(-D__BUILD_GADGETS__)
ENDIF (WIN32)

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
