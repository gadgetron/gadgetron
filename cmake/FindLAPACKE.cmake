
find_package(PkgConfig)
pkg_check_modules(PC_LAPACKE QUIET lapacke)

find_package(LAPACK REQUIRED)

find_path(LAPACKE_INCLUDE_DIR
        NAMES lapacke.h
        PATHS ${PC_LAPACKE_INCLUDE_DIRS}
        )

set(LAPACKE_LIBRARIES ${PC_LAPACKE_LIBRARIES} ${PC_LAPACKE_LINK_LIBRARIES})

message("pkgconfig ${PC_LAPACKE_LIBRARIES} and ${PC_lapacke_LINK_LIBRARIES}")

mark_as_advanced(LAPACKE_FOUND LAPACKE_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACKE
        REQUIRED_VARS LAPACKE_INCLUDE_DIR LAPACKE_LIBRARIES
        )


if(LAPACKE_FOUND)
    add_library(LAPACKE::LAPACKE IMPORTED SHARED)
    set_target_properties(LAPACKE::LAPACKE PROPERTIES
            INCLUDE_DIRECTORIES "${LAPACKE_INCLUDE_DIR}"
            LINK_LIBRARIES "${LAPACKE_LIBRARIES}"
            )
endif()