if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_TOOLBOX_HYPER__)
endif ()


cuda_add_library(gadgetron_toolbox_hyper SHARED 
    CSIOperator.cpp
    CSI_utils.cu
	CSIOperator.h
	CSI_utils.h
	gadgetron_toolbox_hyper_export.h
  )

set_target_properties(gadgetron_toolbox_hyper PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
target_include_directories(gadgetron_toolbox_hyper
		PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(gadgetron_toolbox_hyper 
  gadgetron_toolbox_gpucore 
  gadgetron_toolbox_gpunfft
  ${CUDA_LIBRARIES}
  ${CUDA_CUBLAS_LIBRARIES} 
  )

install(TARGETS gadgetron_toolbox_hyper
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES 
    CSIOperator.h
    CSI_utils.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)