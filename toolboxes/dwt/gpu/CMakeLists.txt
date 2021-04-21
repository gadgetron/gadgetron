if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_GPUDWT__)
endif ()

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif()



cuda_add_library(gadgetron_toolbox_gpudwt SHARED 
    cuNDDWT.h
    cuNDDWT.cu
  )


target_include_directories(gadgetron_toolbox_gpudwt
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(gadgetron_toolbox_gpudwt PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_gpudwt 
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_gpucore 
  ${CUDA_LIBRARIES}
  ${CUDA_CUFFT_LIBRARIES} 
  )



install(TARGETS gadgetron_toolbox_gpudwt
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES
  cuNDDWT.h
	cuDWTOperator.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)