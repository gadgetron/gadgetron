
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was gadgetron-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################


find_package(PythonLibs 3  REQUIRED)
find_package(Boost 1.65.0 COMPONENTS system program_options filesystem timer REQUIRED )
if (Boost_VERSION_STRING VERSION_LESS 1.67.0)
    find_package(Boost 1.65.0 COMPONENTS python3 REQUIRED)
    set(Boost_PYTHON3_TARGET Boost::python3)
    else()
    string(REGEX MATCH "^3\\.([0-9]+)\\.[0-9]+" PYTHON_MINOR_VERSION ${PYTHONLIBS_VERSION_STRING} )
    set(PYTHON_MINOR_VERSION ${CMAKE_MATCH_1})
    find_package(Boost 1.65.0 COMPONENTS "python3${PYTHON_MINOR_VERSION}" REQUIRED)
    set(Boost_PYTHON3_FOUND TRUE)
    set(Boost_PYTHON3_TARGET Boost::python3${PYTHON_MINOR_VERSION})
endif()


add_definitions(-DARMA_DONT_USE_WRAPPER -DARMA_USE_CXX11 -DARMA_64BIT_WORD)
set(GADGETRON_INSTALL_CONFIG_PATH "share/gadgetron/config" )
set(GADGETRON_INSTALL_SCHEMA_PATH "share/gadgetron/schema" )
set(GADGETRON_INSTALL_PYTHON_MODULE_PATH "share/gadgetron/python" )
find_package(ISMRMRD REQUIRED )

if(NOT TARGET gadgetron::gadgetron)
    INCLUDE(${CMAKE_CURRENT_LIST_DIR}/gadgetron-targets.cmake)
endif()
