if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_CPUCORE__)
endif ()


set(header_files ../NDArray.h
                ../complext.h
                ../GadgetronException.h
                ../GadgetronTimer.h
                cpucore_export.h 
                hoNDArray.h
                hoNDArray.hxx
				hoNDArray_iterators.h
                hoNDObjectArray.h
                hoNDArray_utils.h
                hoNDArray_fileio.h
                ho2DArray.h
                ho2DArray.hxx
                ho3DArray.h
                ho3DArray.hxx
                ho4DArray.h
                ho4DArray.hxx
                ho5DArray.h
                ho5DArray.hxx
                ho6DArray.h
                ho6DArray.hxx
                ho7DArray.h
                ho7DArray.hxx 
                hoMatrix.h
                hoMatrix.hxx
                hoNDPoint.h
                hoNDBoundaryHandler.h
                hoNDBoundaryHandler.hxx
                hoNDInterpolator.h
                hoNDInterpolatorNearestNeighbor.hxx
                hoNDInterpolatorLinear.hxx
                hoNDInterpolatorBSpline.hxx )

set(image_files image/hoNDImage.h 
            image/hoNDImage.hxx 
            image/hoNDImageContainer2D.h )

set(algorithm_files algorithm/hoNDBSpline.h algorithm/hoNDBSpline.hxx )

source_group(algorithm FILES ${algorithm_files})
source_group(image FILES ${image_files})

add_library(gadgetron_toolbox_cpucore SHARED
                    hoMatrix.cpp 
                    ${header_files} 
                    ${image_files}  
                    ${algorithm_files} )


target_include_directories(
        gadgetron_toolbox_cpucore
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/toolboxes/core>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/toolboxes/core>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/toolboxes/core/cpu>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/apps/gadgetron>
    $<INSTALL_INTERFACE:include/gadgetron>
)

set_target_properties(gadgetron_toolbox_cpucore  PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_log
  ISMRMRD::ISMRMRD
	boost
  )

install(TARGETS gadgetron_toolbox_cpucore
		EXPORT gadgetron-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
        ${header_files}
        image/hoNDImage.h 
        image/hoNDImage.hxx 
        image/hoNDImageContainer2D.h 
        algorithm/hoNDBSpline.h
        algorithm/hoNDBSpline.hxx 
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

add_subdirectory(math)
add_subdirectory(hostutils)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)