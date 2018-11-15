if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_CPUCORE_MATH__)
endif ()

set(CMAKE_CXX_FLAGS -fpermissive)
include_directories(
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${ARMADILLO_INCLUDE_DIRS}
    ${CBLAS_INCLUDE_DIR})

if (MKL_FOUND)
    include_directories( ${MKL_INCLUDE_DIR} )
    link_directories( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
endif ()

set(cpucore_math_header_files
    cpucore_math_export.h
    hoNDArray_math.h
    hoNDImage_util.h
    hoNDImage_util.hxx
    hoNDArray_linalg.h )

set(cpucore_math_src_files 
    hoNDArray_linalg.cpp )

if (ARMADILLO_FOUND)

    set(cpucore_math_header_files 
        ${cpucore_math_header_files}
        hoNDArray_reductions.h
        hoNDArray_reductions.hxx
        hoArmadillo.h
        hoNDArray_elemwise.h
        cpp_blas.h
         )

    set(cpucore_math_src_files 
        ${cpucore_math_src_files} 
        hoNDArray_reductions.cpp
        hoNDArray_elemwise.cpp
        cpp_blas.cpp)

endif ()
add_library(gadgetron_toolbox_cpucore_math SHARED  ${cpucore_math_src_files} ${cpucore_math_header_files})
set_target_properties(gadgetron_toolbox_cpucore_math PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
if (MKL_FOUND)
    target_link_libraries(gadgetron_toolbox_cpucore_math 
        gadgetron_toolbox_cpucore
        gadgetron_toolbox_log
        ${ARMADILLO_LIBRARIES}
        ${FFTW3_LIBRARIES} 
        ${MKL_LIBRARIES}
        ${CBLAS_LIBRARY}
    )
else ()
    target_link_libraries(gadgetron_toolbox_cpucore_math 
        gadgetron_toolbox_cpucore
        gadgetron_toolbox_log
        ${ARMADILLO_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${FFTW3_LIBRARIES}
        ${CBLAS_LIBRARY}
    )
endif ()

install(TARGETS gadgetron_toolbox_cpucore_math DESTINATION lib COMPONENT main)

install(FILES 
    ${cpucore_math_header_files}
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
