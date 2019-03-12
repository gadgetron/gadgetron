
set(BLA_VENDOR OpenBLAS)
find_package(BLAS)
if (BLAS_FOUND)
    set(CBLAS_LIBRARIES ${BLAS_LIBRARIES})
    find_path(CBLAS_INCLUDE_DIR openblas/cblas.h
        HINTS /usr/include /usr/local/include )
    set(CBLAS_INCLUDE_DIR ${CBLAS_INCLUDE_DIR}/openblas)
else()
    unset(BLA_VENDOR)
    find_package(BLAS REQUIRED)
    find_library(CBLAS_LIBRARIES cblas libcblas
        HINTS /usr/lib/ /usr/lib64)
    find_path(CBLAS_INCLUDE_DIR cblas.h
        HINTS /usr/include /usr/local/include )
endif()

find_package_handle_standard_args(CBLAS DEFAULT_MSG CBLAS_LIBRARIES CBLAS_INCLUDE_DIR)
add_library(BLAS INTERFACE IMPORTED)
set_property(TARGET BLAS PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${BLAS_INCLUDE_DIR} ${CBLAS_INCLUDE_DIR})
set_property(TARGET BLAS PROPERTY INTERFACE_LINK_LIBRARIES ${BLAS_LIBRARIES} ${CBLAS_LIBRARIES})