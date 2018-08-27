find_path(CBLAS_INCLUDE_DIR cblas.h
        HINTS /usr/include /usr/local/include)

find_library(CBLAS_LIBRARY cblas libcblas
        HINTS /usr/lib/ /usr/lib64)

find_package_handle_standard_args(CBLAS DEFAULT_MSG CBLAS_LIBRARY CBLAS_INCLUDE_DIR)
