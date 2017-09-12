if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_bartgadget__)
endif ()

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/klt/cpu 
  ${CMAKE_SOURCE_DIR}/toolboxes/image_io   
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math 
  ${Boost_INCLUDE_DIR}
  )

add_library(gadgetron_bart SHARED 
  bartgadget.h
  bartgadget.cpp
  BART_Recon.xml 
  BART_Recon_cloud.xml
)

set_target_properties(gadgetron_bart PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_bart gadgetron_gadgetbase
  gadgetron_mricore 
  gadgetron_toolbox_mri_core 
  gadgetron_toolbox_log
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_cpufft
  gadgetron_toolbox_cpuklt  
  ${ISMRMRD_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}   
  ${Boost_LIBRARIES}
  )
  
if(ARMADILLO_FOUND)
    target_link_libraries(gadgetron_bart gadgetron_toolbox_cpucore_math )
endif()

install (FILES  bartgadget.h
                DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install (TARGETS gadgetron_bart DESTINATION lib COMPONENT main)

set(GADGETRON_INSTALL_BART_PATH share/gadgetron/bart)
install (FILES Sample_GRAPPA_Recon.txt DESTINATION ${GADGETRON_INSTALL_BART_PATH})

install(FILES BART_Recon.xml BART_Recon_cloud.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
