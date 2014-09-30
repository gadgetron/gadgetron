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

add_library(gadgettools SHARED
  gadgettools_export.h 
  GadgetronSlotContainer.h 
  GadgetronConnector.h 
  GadgetronConnector.cpp 
  GadgetCloudController.h 
  GadgetronCloudConnector.h )

target_link_libraries(gadgettools optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} ${Boost_LIBRARIES})

install(TARGETS gadgettools DESTINATION lib COMPONENT main)

install (FILES 
            GadgetCloudController.h 
            GadgetronCloudConnector.h 
            GadgetronConnector.h
            gadgettools_export.h
            GadgetronSlotContainer.h
            DESTINATION include COMPONENT main)

if(ISMRMRD_FOUND)
  add_subdirectory(ismrmrd)
endif(ISMRMRD_FOUND)
