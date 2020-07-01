# - Find the PLplot libraries
# This module defines
#  PLPLOT_INCLUDE_DIR, the directory for the PLplot headers
#  PLPLOT_LIB_DIR, the directory for the PLplot library files
#  PLPLOT_LIBRARIES, the libraries needed to use PLplot.
#  PLPLOT_FOUND, If false, do not try to use PLplot; if true, the macro definition USE_PLPLOT is added.

if (WIN32)
    if (NOT DEFINED ENV{PLPLOT_PATH})
        set(PLPLOT_PATH "C:/Program Files" CACHE PATH "Where the PLplot are stored")
    endif ()
else ()
    set(PLPLOT_PATH "/usr/include" CACHE PATH "Where the PLplot are stored")
endif ()

if (EXISTS ${PLPLOT_PATH}/plplot)
    set(PLPLOT_FOUND TRUE)
    message("PLplot is found at ${PLPLOT_PATH}/plplot")
else ()
    set(PLPLOT_FOUND FALSE)
    message("PLplot is NOT found ... ")
endif ()

if (PLPLOT_FOUND)
    if (WIN32)
        set(PLPLOT_INCLUDE_DIR "${PLPLOT_PATH}/plplot")
        set(PLPLOT_LIB_DIR "${PLPLOT_PATH}/../lib")
        set(PLPLOT_LIBRARIES plplot)
        set(PLPLOT_LIBRARIES ${PLPLOT_LIBRARIES} plplotcxx)
    else ()
        set(PLPLOT_INCLUDE_DIR "${PLPLOT_PATH}/plplot")
        set(PLPLOT_LIB_DIR "${PLPLOT_PATH}/../lib")

        find_library(PLPLOT_LIB
                NAMES plplot plplotd
                HINTS /usr/lib /usr/lib64 /usr/x86_64-linux-gnu)

        find_library(PLPLOT_CXX_LIB
                NAMES plplotcxx plplotcxxd
                HINTS /usr/lib /usr/lib64 /usr/x86_64-linux-gnu)

        set(PLPLOT_LIBRARIES ${PLPLOT_LIB} ${PLPLOT_CXX_LIB})

        message("PLplot libraries: ${PLPLOT_LIBRARIES}")
    endif ()

    add_definitions(-DUSE_PLPLOT)

endif ()
