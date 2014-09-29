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
 )


INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image )

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

add_library(gadgetron_toolbox_gadgettools SHARED
  gadgettools_export.h 
  GadgetronSlotContainer.h 
  GadgetronConnector.h 
  GadgetronConnector.cpp 
  GadgetCloudController.h 
  GadgetronCloudConnector.h )

set_target_properties(gadgetron_toolbox_gadgettools PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_gadgettools optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} ${Boost_LIBRARIES})

install(TARGETS gadgetron_toolbox_gadgettools DESTINATION lib)

install (FILES 
            GadgetCloudController.h 
            GadgetronCloudConnector.h 
            GadgetronConnector.h
            gadgettools_export.h
            GadgetronSlotContainer.h
            DESTINATION include)

if(ISMRMRD_FOUND)
  add_subdirectory(ismrmrd)
endif(ISMRMRD_FOUND)
