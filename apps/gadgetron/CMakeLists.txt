configure_file(gadgetron_config.in gadgetron_config.h)

if (WIN32)
  add_definitions(-DCROW_MSVC_WORKAROUND)
endif()

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_INSTALL_PREFIX}/include
  ${CMAKE_SOURCE_DIR}/apps/gadgetron
  ${CMAKE_SOURCE_DIR}/toolboxes/cloudbus
  ${CMAKE_SOURCE_DIR}/toolboxes/rest
  ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
  ${Boost_INCLUDE_DIR}
  ${ACE_INCLUDE_DIR}
  )

if (CUDA_FOUND)
  include_directories(${CUDA_INCLUDE_DIRS})
endif()

add_executable(gadgetron 
  main.cpp 
  GadgetServerAcceptor.h
  GadgetServerAcceptor.cpp 
  GadgetStreamController.h
  # EndGadget.h
  # EndGadget.cpp
  # GadgetStreamInterface.cpp
  Gadget.h 
  GadgetContainerMessage.h 
  GadgetMessageInterface.h
  GadgetronExport.h 
  gadgetron_xml.h
  )

target_link_libraries(gadgetron 
  gadgetron_gadgetbase
  gadgetron_toolbox_log
  gadgetron_toolbox_rest
  gadgetron_toolbox_gadgettools gadgetron_toolbox_cloudbus 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
 )

if (BUILD_PYTHON_SUPPORT)
    add_definitions("-DCOMPILING_WITH_PYTHON_SUPPORT")
endif()

add_executable(gadgetron_info
  gadgetron_info.cpp
  EndGadget.cpp
  GadgetStreamInterface.cpp
)

target_link_libraries(gadgetron_info 
  gadgetron_toolbox_log
  gadgetron_gadgetbase
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
 )

if (CUDA_FOUND)
  target_link_libraries(gadgetron_info ${CUDA_LIBRARIES})
  target_link_libraries(gadgetron ${CUDA_LIBRARIES})
endif()

add_library(gadgetron_gadgetbase SHARED
  Gadget.h
  Gadget.cpp
  GadgetStreamController.h
  GadgetStreamController.cpp
  gadgetron_xml.h
  gadgetron_xml.cpp
  gadgetron_home.h
  gadgetron_home.cpp
  pugixml.cpp  
  GadgetStreamInterface.h 
  GadgetStreamInterface.cpp 
  EndGadget.h
  EndGadget.cpp
)

target_link_libraries(gadgetron_gadgetbase
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  gadgetron_toolbox_gadgettools
  gadgetron_toolbox_cloudbus
  gadgetron_toolbox_log
)

if (REQUIRE_SIGNED_CONFIG)
    target_link_libraries(gadgetron_info GTBabylon)
    target_link_libraries(gadgetron GTBabylon)
    target_link_libraries(gadgetron_gadgetbase GTBabylon)
endif()

if (GPERFTOOLS_PROFILER)
  message("Adding gperftools cpu profiler to gadgetron link assemblage.")
  target_link_libraries(gadgetron ${GPERFTOOLS_PROFILER} ${GPERFTOOLS_TCMALLOC})
endif()


set_target_properties(gadgetron_gadgetbase PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
set_target_properties (gadgetron_gadgetbase PROPERTIES COMPILE_DEFINITIONS "__BUILD_GADGETRON_GADGETBASE__")

install(TARGETS gadgetron gadgetron_info DESTINATION bin COMPONENT main)
install(TARGETS gadgetron_gadgetbase DESTINATION lib COMPONENT main)

install(FILES
  gadgetbase_export.h
  EndGadget.h
  Gadget.h
  GadgetContainerMessage.h
  GadgetMessageInterface.h
  GadgetronExport.h
  gadgetron_home.h
  gadgetron_xml.h
  GadgetServerAcceptor.h
  GadgetStreamController.h
  GadgetStreamInterface.h
  gadgetron_home.h
  ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_config.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main) 

install(FILES 
  gadgetron.xml.example
  DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

install(FILES 
  schema/gadgetron.xsd
  DESTINATION ${GADGETRON_INSTALL_SCHEMA_PATH} COMPONENT main)

