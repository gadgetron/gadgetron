
IF (WIN32)
    ADD_DEFINITIONS(-DTIXML_USE_STL)
ENDIF (WIN32)

find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

include_directories(${ACE_INCLUDE_DIR} 
                    ${Boost_INCLUDE_DIR} 
                    ${XSD_INCLUDE_DIR} 
                    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
                    ${CMAKE_SOURCE_DIR}/toolboxes/core
                    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
                    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
                    ${CMAKE_SOURCE_DIR}/apps/gadgetron
                    ${CMAKE_SOURCE_DIR}/gadgets/mri_core 
                    ${ISMRMRD_XSD_INCLUDE_DIR} 
                    ${ISMRMRD_INCLUDE_DIR} )

#Process the XSD files
SET(XSDS schema/gadgetron.xsd)
SET(XSD_ARGS cxx-tree --generate-serialization)
WRAP_XSD(XSDS_SOURCES XSD_INCLUDES ${CMAKE_CURRENT_BINARY_DIR}/schema ${XSDS} OPTIONS ${XSD_ARGS})
INCLUDE_DIRECTORIES(${XSD_INCLUDES} ${XERCESC_INCLUDE_DIR})

message("XSDS_SOURCES is " ${XSDS_SOURCES})

if (MKL_FOUND)
    INCLUDE_DIRECTORIES( ${MKL_INCLUDE_DIR} )
endif (MKL_FOUND)

add_library(gadgettools SHARED ${XSDS_SOURCES} 
                               GadgetImageMessageReader.h 
                               GadgetImageMessageWriter.h 
                               gadgettools_export.h 
                               GadgetronSlotContainer.h 
                               GadgetronConnector.h 
                               GadgetronConnector.cpp 
                               GadgetServerAcceptor.h 
                               GadgetServerAcceptor.cpp 
                               GadgetStreamController.h 
                               GadgetStreamController.cpp
                               GadgetCloudController.h 
                               GadgetronCloudConnector.h )

target_link_libraries(gadgettools optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} ${XERCESC_LIBRARIES})

install(TARGETS gadgettools DESTINATION lib)

install (FILES  GadgetImageMessageReader.h 
            GadgetImageMessageWriter.h
            GadgetronConnector.h
            gadgettools_export.h
            GadgetronSlotContainer.h
            GadgetServerAcceptor.h
            GadgetStreamController.h
            DESTINATION include)
