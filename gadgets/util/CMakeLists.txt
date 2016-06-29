if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_UTIL_GADGETS__)
endif (WIN32)

include_directories(
    ${CMAKE_SOURCE_DIR}/toolboxes/core
)

add_library(gadgetron_util SHARED 
    gadgetron_util_gadgets_export.h 
    ParameterRelayGadget.h
    ParameterRelayGadget.cpp
    GadgetMessageTextWriter.h
    GadgetMessageTextWriter.cpp
    GadgetronInfoGadget.h
    GadgetronInfoGadget.cpp
)

set_target_properties(gadgetron_util PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})                                                                                                                                                                                                      
if (CUDA_FOUND)
target_link_libraries(gadgetron_util
    gadgetron_gadgetbase
    gadgetron_toolbox_log
    ${ACE_LIBRARIES}
    ${CUDA_LIBRARIES} 
    ${CUDA_CUFFT_LIBRARIES} 
    ${CUDA_CUBLAS_LIBRARIES} 
)

else (CUDA_FOUND)
  target_link_libraries(gadgetron_util
    gadgetron_gadgetbase
    gadgetron_toolbox_log
    ${ACE_LIBRARIES}
)
endif(CUDA_FOUND)

install(FILES 
    gadgetron_util_gadgets_export.h
    ParameterRelayGadget.h
    GadgetMessageTextWriter.h
    GadgetronInfoGadget.h
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES config/gadgetron_info.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
install(TARGETS gadgetron_util DESTINATION lib COMPONENT main)
