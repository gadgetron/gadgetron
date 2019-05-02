find_path(CBLAS_INCLUDE_DIR cblas.h
        HINTS /usr/include /usr/local/include /usr/include/openblas /usr/local/include/openblas)

find_library(CBLAS_LIBRARY cblas libcblas openblas libopenblas
        HINTS /usr/lib/ /usr/lib64)

find_package_handle_standard_args(CBLAS DEFAULT_MSG CBLAS_LIBRARY CBLAS_INCLUDE_DIR)