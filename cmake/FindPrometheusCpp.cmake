#
# Find the prometheus-cpp library
#    https://github.com/jupp0r/prometheus-cpp
#

# This module defines
# PROMETHEUSCPP_INCLUDE_DIR, where to find exposer.h
# PROMETHEUSCPP_LIBRARIES, the libraries to link against
# PROMETHEUSCPP_FOUND, if false, you cannot build anything that requires ZFP

######################################################################## 
find_path(PROMETHEUSCPP_INCLUDE_DIR prometheus/exposer.h /usr/include /usr/local/include  DOC "directory containing promeths/exposer.h for PROMETHEUSCPP library")  
find_library(PROMETHEUS_CORE_LIBRARY NAMES prometheus-cpp-core PATHS /usr/lib /usr/local/lib DOC "PROMETHEUS CORE library file") 
find_library(PROMETHEUS_PULL_LIBRARY NAMES prometheus-cpp-pull PATHS /usr/lib /usr/local/lib DOC "PROMETHEUS PULL library file") 
find_library(PROMETHEUS_PUSH_LIBRARY NAMES prometheus-cpp-push PATHS /usr/lib /usr/local/lib DOC "PROMETHEUS PUSH library file") 
 
if (WIN32) 
    set(PROMETHEUS_CPP_FOUND FALSE)
else ()
    if(PROMETHEUSCPP_INCLUDE_DIR AND PROMETHEUS_CORE_LIBRARY)
        set(PROMETHEUSCPP_LIBRARIES ${PROMETHEUS_CORE_LIBRARY} ${PROMETHEUS_PULL_LIBRARY} ${PROMETHEUS_PUSH_LIBRARY})
        set(PROMETHEUSCPP_FOUND TRUE)
    else ()
        set(PROMETHEUSCPP_FOUND FALSE)
    endif ()
endif () 
 
if (PROMETHEUSCPP_FOUND)
    message("PrometheusCpp Found")
    message("PROMETHEUSCPP_INCLUDE_DIR ${PROMETHEUSCPP_INCLUDE_DIR}")
    message("PROMETHEUS_LIBRARIES ${PROMETHEUSCPP_LIBRARIES}")
else ()
    message("PrometheusCpp not found")
endif ()