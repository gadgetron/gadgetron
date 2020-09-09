if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_GPUSDC__)
  add_definitions(-D_USE_MATH_DEFINES)
endif ()

cuda_add_library(gadgetron_toolbox_gpusdc SHARED 
  cuSDC.cu
)

set_target_properties(gadgetron_toolbox_gpusdc PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_include_directories(gadgetron_toolbox_gpusdc
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(gadgetron_toolbox_gpusdc
  gadgetron_toolbox_gpucore
  gadgetron_toolbox_operator
  gadgetron_toolbox_log
  gadgetron_toolbox_sdc
  gadgetron_toolbox_nfft
  gadgetron_toolbox_gpunfft
  ${CUDA_LIBRARIES} 
)

install(TARGETS gadgetron_toolbox_gpusdc
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES 
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH}
  COMPONENT main)
