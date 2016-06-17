if (WIN32)
    ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUCORE__)
endif (WIN32)

include_directories(
    ${FFTW3_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/apps/gadgetron
)

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

#if (MKL_FOUND)
#    INCLUDE_DIRECTORIES( ${MKL_INCLUDE_DIR} )
#    LINK_DIRECTORIES( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
#endif (MKL_FOUND)

set(header_files ../NDArray.h
                ../complext.h
                ../GadgetronException.h
                ../GadgetronTimer.h
                cpucore_export.h 
                hoNDArray.h
                hoNDArray.hxx
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

set(algorithm_files algorithm/hoNDBSpline.h algorithm/hoNDBSpline.hxx)

source_group(algorithm FILES ${algorithm_files})
source_group(image FILES ${image_files})

add_library(gadgetron_toolbox_cpucore SHARED
                    hoMatrix.cpp 
                    ${header_files} 
                    ${image_files}  
                    ${algorithm_files} )

set_target_properties(gadgetron_toolbox_cpucore  PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_log
  )

install(TARGETS gadgetron_toolbox_cpucore DESTINATION lib COMPONENT main)

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
