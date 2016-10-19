if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_BART_Gadget__)
endif ()

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  )

add_library(gadgetron_BART_Gadget SHARED 
  BART_Gadget.h
  BART_Gadget.cpp
)

set_target_properties(gadgetron_BART_Gadget PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_BART_Gadget gadgetron_toolbox_cpucore 
  gadgetron_gadgetbase
  gadgetron_toolbox_log
  gadgetron_toolbox_cpufft
  ${ISMRMRD_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

find_package(Boost COMPONENTS tokenizer) 
if(Boost_FOUND)
    target_link_libraries(gadgetron_BART_Gadget ${Boost_LIBRARIES})
endif()

if(ARMADILLO_FOUND)
    target_link_libraries(gadgetron_BART_Gadget gadgetron_toolbox_cpucore_math )
endif()

install (FILES  BART_Gadget.h
                DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install (TARGETS gadgetron_BART_Gadget DESTINATION lib COMPONENT main)

install (FILES BART_Commands__.txt DESTINATION lib COMPONENT main)

install(FILES
  BART_Recon.xml
  DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
