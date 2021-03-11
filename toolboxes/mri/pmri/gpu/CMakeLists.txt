if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_GPUPMRI__)
  add_definitions(-DWIN32_LEAN_AND_MEAN)
endif ()

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif()


cuda_add_library(gadgetron_toolbox_gpuparallelmri SHARED 
    b1_map.h
    spirit_calibration.h
    cuCartesianSenseOperator.h
    cuNonCartesianKtSenseOperator.h
    cuNonCartesianSenseOperator.h
    cuSpiritOperator.h
    cuBuffer.h
    cuSenseBuffer.h
    cuSenseBufferCg.h
    cuSenseOperator.h
    gpupmri_export.h
    htgrappa.h
    senseOperator.h
    sense_utilities.h
    b1_map.cu
    b1_map_NIH_Souheil.cu
    spirit_calibration.cu
    sense_utilities.cu
    cuSenseOperator.cu
    cuCartesianSenseOperator.cu
    cuNonCartesianSenseOperator.cu
    cuNonCartesianKtSenseOperator.cu
    cuBuffer.cpp
    cuSenseBuffer.cpp
    cuSenseBufferCg.cpp
    cuSpiritBuffer.cpp
    htgrappa.cpp
    htgrappa.cu
    trajectory_utils.h
    trajectory_utils.cu
  )

set_target_properties(gadgetron_toolbox_gpuparallelmri PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_gpuparallelmri 
  gadgetron_toolbox_gpucore
  gadgetron_toolbox_log
  gadgetron_toolbox_gpusolvers
  gadgetron_toolbox_gpunfft 
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_cpucore_math
  )


target_include_directories(gadgetron_toolbox_gpuparallelmri
		PUBLIC
		$<INSTALL_INTERFACE:include>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
		)

install(TARGETS gadgetron_toolbox_gpuparallelmri
    EXPORT gadgetron-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
	b1_map.h
	sense_utilities.h
	htgrappa.h
	senseOperator.h
	cuSenseOperator.h
	cuCartesianSenseOperator.h
	cuNonCartesianSenseOperator.h
	cuNonCartesianKtSenseOperator.h
        cuSpiritOperator.h
        cuBuffer.h
	cuSenseBuffer.h
	cuSenseBufferCg.h
	cuSpiritBuffer.h
	gpupmri_export.h
DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)