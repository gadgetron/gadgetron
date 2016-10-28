if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_bartgadget__)
endif ()

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  )

add_library(gadgetron_bart SHARED 
  bartgadget.h
  bartgadget.cpp
)

set_target_properties(gadgetron_bart PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_bart gadgetron_toolbox_cpucore 
  gadgetron_gadgetbase
  gadgetron_toolbox_log
  gadgetron_toolbox_cpufft
  ${ISMRMRD_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

find_package(Boost COMPONENTS filesystem) 
if(Boost_FOUND)
    target_link_libraries(gadgetron_bart ${Boost_LIBRARIES})
endif()

if(ARMADILLO_FOUND)
    target_link_libraries(gadgetron_bart gadgetron_toolbox_cpucore_math )
endif()

install (FILES  bartgadget.h
                DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install (TARGETS gadgetron_bart DESTINATION lib COMPONENT main)

set(GADGETRON_INSTALL_BART_PATH share/gadgetron/bart)
install (FILES Sample_GRAPPA_Recon.txt DESTINATION ${GADGETRON_INSTALL_BART_PATH})

install(FILES BART_Recon.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)