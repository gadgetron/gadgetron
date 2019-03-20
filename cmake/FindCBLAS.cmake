if (VCPKG_TARGET_TRIPLET)
find_package(OpenBLAS CONFIG REQUIRED)
add_library(BLAS INTERFACE IMPORTED)

set_property(TARGET BLAS PROPERTY INTERFACE_LINK_LIBRARIES OpenBLAS::OpenBLAS)
else()

set(BLA_VENDOR OpenBLAS)
find_package(BLAS)
if (BLAS_FOUND)
message("OpenBLAS found")
    set(CBLAS_LIBRARIES ${BLAS_LIBRARIES})
    find_path(CBLAS_INCLUDE_DIR openblas/cblas.h
        PATHS /usr/include /usr/local/include )
    set(CBLAS_INCLUDE_DIR ${CBLAS_INCLUDE_DIR}/openblas)
else()
    unset(BLA_VENDOR)
    find_package(BLAS REQUIRED)
    find_library(CBLAS_LIBRARIES cblas libcblas
        PATHS /usr/lib/ /usr/lib64)
    find_path(CBLAS_INCLUDE_DIR cblas.h
        PATHS /usr/include /usr/local/include )
endif()

find_package_handle_standard_args(CBLAS DEFAULT_MSG CBLAS_LIBRARIES CBLAS_INCLUDE_DIR)
add_library(BLAS INTERFACE IMPORTED)
set_property(TARGET BLAS PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${BLAS_INCLUDE_DIR} ${CBLAS_INCLUDE_DIR})
set_property(TARGET BLAS PROPERTY INTERFACE_LINK_LIBRARIES ${BLAS_LIBRARIES} ${CBLAS_LIBRARIES})
endif()