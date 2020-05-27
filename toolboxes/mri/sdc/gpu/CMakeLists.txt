if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_GPUSDC__)
  add_definitions(-D_USE_MATH_DEFINES)
endif ()

cuda_add_library(gadgetron_toolbox_gpusdc SHARED 
  cuSDC.h
  cuSDC.cu
  cuSDC_kernel.h
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
  ${CUDA_LIBRARIES} 
)

install(TARGETS gadgetron_toolbox_gpusdc
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES 
  cuSDC.h 
  cuSDC_kernel.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH}
  COMPONENT main)
