if (WIN32)
    ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUCORE_MATH__)
endif (WIN32)

include_directories(
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${ARMADILLO_INCLUDE_DIRS} 
    ${ISMRMRD_INCLUDE_DIR} )

if (MKL_FOUND)
    INCLUDE_DIRECTORIES( ${MKL_INCLUDE_DIR} )
    LINK_DIRECTORIES( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
endif (MKL_FOUND)

set(cpucore_math_header_files
    cpucore_math_export.h
    # hoNDArray_math_util.h
    hoNDArray_math.h
    hoNDImage_util.h
    hoNDImage_util.hxx
    # hoNDImage_util_instantiate.hxx 
    hoNDArray_linalg.h )

set(cpucore_math_src_files 
    # hoNDArray_math_util.cpp
    # hoNDImage_util.cpp 
    hoNDArray_linalg.cpp )

if (ARMADILLO_FOUND)

    set(cpucore_math_header_files 
        ${cpucore_math_header_files}
        hoNDArray_reductions.h
        hoArmadillo.h
        hoNDArray_elemwise.h
         )

    set(cpucore_math_src_files 
        ${cpucore_math_src_files} 
        hoNDArray_reductions.cpp
        hoNDArray_elemwise.cpp
        )

endif (ARMADILLO_FOUND)

add_library(gadgetron_toolbox_cpucore_math SHARED ${cpucore_math_header_files} ${cpucore_math_src_files})
set_target_properties(gadgetron_toolbox_cpucore_math PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

if (MKL_FOUND)
    target_link_libraries(gadgetron_toolbox_cpucore_math 
        gadgetron_toolbox_cpucore
	gadgetron_toolbox_log
        ${ARMADILLO_LIBRARIES}
        ${ISMRMRD_LIBRARIES}
        ${FFTW3_LIBRARIES} 
        ${MKL_LIBRARIES}
    )
else (MKL_FOUND)
    target_link_libraries(gadgetron_toolbox_cpucore_math 
        gadgetron_toolbox_cpucore
	gadgetron_toolbox_log
        ${ARMADILLO_LIBRARIES}
        ${ISMRMRD_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${FFTW3_LIBRARIES}
    )
endif (MKL_FOUND)

install(TARGETS gadgetron_toolbox_cpucore_math DESTINATION lib COMPONENT main)

install(FILES 
    ${cpucore_math_header_files}
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
