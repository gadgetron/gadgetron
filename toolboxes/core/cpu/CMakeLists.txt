if (WIN32)
    ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUCORE__)
endif (WIN32)

find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

include_directories(
    ${XSD_INCLUDE_DIR}
    ${FFTW3_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
    ${CMAKE_SOURCE_DIR}/apps/gadgetron
)

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

#if (MKL_FOUND)
#    INCLUDE_DIRECTORIES( ${MKL_INCLUDE_DIR} )
#    LINK_DIRECTORIES( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
#endif (MKL_FOUND)

#Process the XSD files
SET(XSDS image/schema/gadgetronMetaAttributes.xsd)
SET(XSD_ARGS cxx-tree --generate-serialization)
WRAP_XSD(XSDS_SOURCES XSD_INCLUDES ${CMAKE_CURRENT_BINARY_DIR}/schema ${XSDS} OPTIONS ${XSD_ARGS})
INCLUDE_DIRECTORIES(${XSD_INCLUDES} ${XERCESC_INCLUDE_DIR})

message("XSDS_SOURCES is " ${XSDS_SOURCES})

set(header_files ../NDArray.h
                ../complext.h
                ../GadgetronException.h
                ../GadgetronCommon.h
                ../GadgetronTimer.h
                ../SerializableObject.h
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
            image/hoNDImageContainer2D.h 
            image/hoNDMetaAttributesBase.h 
            image/hoNDMetaAttributes.h 
            image/hoNDMetaAttributes.cpp 
            image/schema/gadgetronMetaAttributes.xsd )

set(algorithm_files algorithm/hoNDBSpline.h algorithm/hoNDBSpline.hxx)

source_group(algorithm FILES ${algorithm_files})
source_group(image FILES ${image_files})

add_library(cpucore ${LIBTYPE}
                    ${XSDS_SOURCES} 
                    hoMatrix.cpp 
                    ${header_files} 
                    ${image_files}  
                    ${algorithm_files} )

target_link_libraries(cpucore 
                    ${Boost_LIBRARIES} 
                    ${XERCESC_LIBRARIES} )

install(TARGETS cpucore DESTINATION lib)

install(FILES
        ${header_files}
        image/hoNDImage.h 
        image/hoNDImage.hxx 
        image/hoNDImageContainer2D.h 
        image/hoNDMetaAttributesBase.h 
        image/hoNDMetaAttributes.h 
        algorithm/hoNDBSpline.h
        algorithm/hoNDBSpline.hxx 
        DESTINATION include)

install(FILES image/schema/gadgetronMetaAttributes.xsd DESTINATION schema)

add_subdirectory(arma_math)
add_subdirectory(hostutils)
