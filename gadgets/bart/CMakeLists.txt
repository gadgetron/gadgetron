
cmake_policy(SET CMP0028 NEW) # double colon only for imported or alias libraries

# ==============================================================================
# First try to find either the BARTmain library or a path to the BART source code

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

find_package(BART REQUIRED CONFIG)

# ==============================================================================

set(GADGET_FILES
  bartgadget.h
  bartgadget.cpp
  bart_helpers.h
  bart_helpers.cpp
  bart_logger.cpp
  BART_Recon.xml
  BART_Recon_cloud.xml
  BART_Recon_cloud_Standard.xml
)

# ------------------------------------------------------------------------------

if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_bartgadget__)
endif ()

if(BART_MEMONLY_CFL)
  add_definitions(-DMEMONLY_CFL)
endif()

add_library(gadgetron_bart SHARED ${GADGET_FILES})
target_link_libraries(gadgetron_bart BART::BARTMAIN)



set_target_properties(gadgetron_bart
  PROPERTIES
  VERSION ${GADGETRON_VERSION_STRING}
  SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_bart gadgetron_core
  gadgetron_mricore 
  gadgetron_toolbox_mri_core 
  gadgetron_toolbox_log
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_cpufft
  gadgetron_toolbox_cpuklt  
  )

  target_link_libraries(gadgetron_bart gadgetron_toolbox_cpucore_math)

install(FILES  bartgadget.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

  install(TARGETS gadgetron_bart
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

set(GADGETRON_INSTALL_BART_PATH share/gadgetron/bart)
install(FILES Sample_Grappa_Recon.sh Sample_Grappa_Recon_Standard.sh DESTINATION ${GADGETRON_INSTALL_BART_PATH})

install(FILES BART_Recon.xml BART_Recon_cloud.xml BART_Recon_cloud_Standard.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)