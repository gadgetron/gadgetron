# - Find the PLplot libraries
# This module defines
#  PLPLOT_INCLUDE_DIR, the directory for the PLplot headers
#  PLPLOT_LIB_DIR, the directory for the PLplot library files
#  PLPLOT_LIBRARIES, the libraries needed to use PLplot.
#  PLPLOT_FOUND, If false, do not try to use PLplot; if true, the macro definition USE_PLPLOT is added.

if ( WIN32 )
  if(NOT DEFINED ENV{PLPLOT_PATH})
    set(PLPLOT_PATH "C:/Program Files" CACHE PATH "Where the PLplot are stored")
  endif() 
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
    if ( WIN32 )
        set(PLPLOT_INCLUDE_DIR "${PLPLOT_PATH}/plplot/include/plplot")
        set(PLPLOT_LIB_DIR "${PLPLOT_PATH}/plplot/lib")
        set(PLPLOT_LIBRARIES plplot)
        set(PLPLOT_LIBRARIES ${PLPLOT_LIBRARIES} plplotcxx)
    else ()
        set(PLPLOT_INCLUDE_DIR "${PLPLOT_PATH}/plplot")
        set(PLPLOT_LIB_DIR "${PLPLOT_PATH}/../lib")
        # linux distribution has older PLplot version
        # the 'd' means double precision, not debug
        set(PLPLOT_LIBRARIES plplotd)
        set(PLPLOT_LIBRARIES ${PLPLOT_LIBRARIES} plplotcxxd)
    endif ()

    add_definitions(-DUSE_PLPLOT)

endif ()
