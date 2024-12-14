# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
# TODO: DELETE ^^^^^

set(cpucore_math_header_files
    hoNDArray_math.h
    hoNDImage_util.h
    hoNDImage_util.hxx
    hoNDArray_linalg.h
    hoNDArray_reductions.h
    hoNDArray_reductions.hxx
    hoArmadillo.h
    hoNDArray_elemwise.h
    hoNDArray_elemwise.hpp
    cpp_blas.h
    cpp_lapack.h
    )

set(cpucore_math_src_files
    hoNDArray_linalg.cpp
    hoNDArray_reductions.cpp
    hoNDArray_elemwise.cpp
    cpp_blas.cpp
    cpp_lapack.cpp
    )

set_source_files_properties(cpp_blas.cpp PROPERTIES COMPILE_FLAGS -fpermissive)
add_library(pingvin_toolbox_cpucore_math SHARED  ${cpucore_math_src_files} ${cpucore_math_header_files})
set_target_properties(pingvin_toolbox_cpucore_math PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})
target_link_libraries(pingvin_toolbox_cpucore_math
    pingvin_toolbox_cpucore
    pingvin_toolbox_log
    armadillo
    GTBLAS
    ${LAPACK_LIBRARIES}
        )

install(TARGETS pingvin_toolbox_cpucore_math
        EXPORT pingvin-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES
    ${cpucore_math_header_files}
    DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)