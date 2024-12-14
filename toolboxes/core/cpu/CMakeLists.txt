set(header_files
                hoNDArray.h
                hoNDArray.hxx
                hoNDArray_converter.h
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

add_library(pingvin_toolbox_cpucore SHARED
                hoMatrix.cpp
                ../NDArray.h
                ../complext.h
                ../GadgetronTimer.h
                ${header_files}
                ${image_files}
                ${algorithm_files}
)

target_include_directories(
        pingvin_toolbox_cpucore
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/core>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/toolboxes/core>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/toolboxes/core>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/toolboxes/core/cpu>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math>
    $<INSTALL_INTERFACE:${PINGVIN_INSTALL_INCLUDE_PATH}>
)

set_target_properties(pingvin_toolbox_cpucore  PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(
  pingvin_toolbox_cpucore
  pingvin_toolbox_log
	Boost::boost
    Boost::system
    Boost::filesystem
  )

install(TARGETS pingvin_toolbox_cpucore
		EXPORT pingvin-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
        ${header_files}
        ${image_files}
        ${algorithm_files}
        DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)

add_subdirectory(math)