
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

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/klt/cpu 
  ${CMAKE_SOURCE_DIR}/toolboxes/image_io   
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math 
  ${Boost_INCLUDE_DIR}
  )

add_library(gadgetron_bart SHARED ${GADGET_FILES})
target_link_libraries(gadgetron_bart BART::BARTMAIN)

if(USE_CUDA)
  find_package(CUDA)
  CUDA_ADD_CUFFT_TO_TARGET(gadgetron_bart)
  CUDA_ADD_CUBLAS_TO_TARGET(gadgetron_bart)
  target_link_libraries(gadgetron_bart ${CUDA_LIBRARIES})
endif()


set_target_properties(gadgetron_bart
  PROPERTIES
  VERSION ${GADGETRON_VERSION_STRING}
  SOVERSION ${GADGETRON_SOVERSION})

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
  target_link_libraries(gadgetron_bart gadgetron_toolbox_cpucore_math)
endif()

install(FILES  bartgadget.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_bart DESTINATION lib COMPONENT main)

set(GADGETRON_INSTALL_BART_PATH share/gadgetron/bart)
install(FILES Sample_Grappa_Recon.sh Sample_Grappa_Recon_Standard.sh DESTINATION ${GADGETRON_INSTALL_BART_PATH})

install(FILES BART_Recon.xml BART_Recon_cloud.xml BART_Recon_cloud_Standard.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
