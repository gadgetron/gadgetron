find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

include_directories(
                    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
                    ${CMAKE_SOURCE_DIR}/toolboxes/core
                    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
                    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
                    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
                    ${CMAKE_SOURCE_DIR}/apps/gadgetron
                    ${CMAKE_SOURCE_DIR}/gadgets/mri_core 
                    ${ACE_INCLUDE_DIR} 
                    ${Boost_INCLUDE_DIR} 
                    ${XSD_INCLUDE_DIR} 
 )

#Process the XSD files
SET(XSDS schema/gadgetron.xsd)
SET(XSD_ARGS cxx-tree --generate-serialization)
WRAP_XSD(XSDS_SOURCES XSD_INCLUDES ${CMAKE_CURRENT_BINARY_DIR}/schema ${XSDS} OPTIONS ${XSD_ARGS})
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
                    ${XSD_INCLUDES} 
                    ${XERCESC_INCLUDE_DIR} )

message("XSDS_SOURCES is " ${XSDS_SOURCES})

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

add_library(gadgettools SHARED ${XSDS_SOURCES} 
                               gadgettools_export.h 
                               GadgetronSlotContainer.h 
                               GadgetronConnector.h 
                               GadgetronConnector.cpp 
                               GadgetCloudController.h 
                               GadgetronCloudConnector.h 
                               schema/gadgetron.xsd )

target_link_libraries(gadgettools optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} ${XERCESC_LIBRARIES} ${Boost_LIBRARIES})

install(TARGETS gadgettools DESTINATION lib)

install (FILES 
            GadgetronConnector.h
            gadgettools_export.h
            GadgetronSlotContainer.h
            DESTINATION include)

if(ISMRMRD_FOUND)
  add_subdirectory(ismrmrd)
endif(ISMRMRD_FOUND)
