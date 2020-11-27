if (VCPKG_TARGET_TRIPLET)
    find_package(OpenBLAS CONFIG REQUIRED)
    add_library(BLAS INTERFACE IMPORTED)

    set_property(TARGET BLAS PROPERTY INTERFACE_LINK_LIBRARIES OpenBLAS::OpenBLAS)
	set(BLAS_LIBRARIES OpenBLAS::OpenBLAS)
else ()

    set(BLA_VENDOR OpenBLAS)
    find_package(BLAS)
    if (BLAS_FOUND)
        message("OpenBLAS found")

        ##Let's look for OpenBLAS OpenMP
        find_library(OPENBLAS_OMP_LIBRARIES openblaso)
        message("OPENBLAS_OMP_LIBRARIES ${OPENBLAS_OMP_LIBRARIES} ${OPENBLAS_OMP_LIBRARIES_FOUND}")
        if (NOT OPENBLAS_OMP_LIBRARIES)
            message("OMP OPENBLAS NOT FOUND")
            set(CBLAS_LIBRARIES ${BLAS_LIBRARIES})
        else()
            message("OMP OPENBLAS FOUND")
            set(CBLAS_LIBRARIES ${OPENBLAS_OMP_LIBRARIES})
        endif()
        find_path(CBLAS_BASE_INCLUDE_DIR openblas/cblas.h
                PATHS /usr/include /usr/local/include)
        if (NOT CBLAS_BASE_INCLUDE_DIR)
            find_path(CBLAS_INCLUDE_DIR cblas.h
                    PATHS /usr/include /usr/local/include)
        else ()
            set(CBLAS_INCLUDE_DIR ${CBLAS_BASE_INCLUDE_DIR}/openblas)
        endif ()
    else ()
        unset(BLA_VENDOR)
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)

        set(CBLAS_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
        find_path(CBLAS_INCLUDE_DIR cblas.h
                  PATHS /usr/include /usr/local/include)
    endif ()

    message("BLAS ${CBLAS_LIBRARIES} LAPACK ${LAPACK_LIBRARIES}" )

    find_package_handle_standard_args(CBLAS FOUND_VAR CBLAS_FOUND REQUIRED_VARS CBLAS_LIBRARIES CBLAS_INCLUDE_DIR)

endif ()