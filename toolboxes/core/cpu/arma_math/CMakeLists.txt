if (WIN32)
    ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUCORE_MATH__)
endif (WIN32)

include_directories(
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
    ${ARMADILLO_INCLUDE_DIRS} )

if (MKL_FOUND)
    INCLUDE_DIRECTORIES( ${MKL_INCLUDE_DIR} )
    LINK_DIRECTORIES( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
endif (MKL_FOUND)

set(cpucore_math_header_files
    cpucore_math_export.h
    hoNDArray_math_util.h
    hoNDArray_math_util.hxx
    hoNDImage_util.h
    hoNDImage_util.hxx
    hoNDImage_util_instantiate.hxx 
    hoMatrix_util.h )

set(cpucore_math_src_files 
    hoNDArray_math_util.cpp
    hoNDImage_util.cpp 
    hoMatrix_util.cpp )

if (ARMADILLO_FOUND)

    set(cpucore_math_header_files 
        ${cpucore_math_header_files}
        hoNDArray_operators.h
        hoNDArray_reductions.h
        hoArmadillo.h
        hoNDArray_elemwise.h
        hoNDArray_blas.h )

    set(cpucore_math_src_files 
        ${cpucore_math_src_files} 
        hoNDArray_operators.cpp
        hoNDArray_reductions.cpp
        hoNDArray_elemwise.cpp
        hoNDArray_blas.cpp )

endif (ARMADILLO_FOUND)

add_library(cpucore_math ${LIBTYPE} ${cpucore_math_header_files} ${cpucore_math_src_files})

message(STATUS "arma: ${ARMADILLO_LIBRARIES}")

if (MKL_FOUND)
    target_link_libraries(cpucore_math 
        cpucore
        ${ARMADILLO_LIBRARIES}
        ${FFTW3_LIBRARIES} 
        ${MKL_LIBRARIES}
    )
else (MKL_FOUND)
    target_link_libraries(cpucore_math 
        cpucore
        ${ARMADILLO_LIBRARIES}
        ${FFTW3_LIBRARIES}
    )
endif (MKL_FOUND)

install(TARGETS cpucore_math DESTINATION lib)

install(FILES 
    ${cpucore_math_header_files}
    DESTINATION include)
