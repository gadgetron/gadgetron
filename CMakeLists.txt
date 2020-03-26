cmake_minimum_required(VERSION 3.10...3.12 FATAL_ERROR)

if (${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else ()
    cmake_policy(VERSION 3.12)
endif ()

project(GADGETRON LANGUAGES CXX C)

find_program(CCACHE_FOUND ccache)
if (CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif (CCACHE_FOUND)

if (EXISTS $ENV{CMAKE_PREFIX_PATH})
    list(APPEND CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH})
endif ()
if (EXISTS $ENV{ISMRMRD_HOME})
    list(APPEND CMAKE_PREFIX_PATH $ENV{ISMRMRD_HOME})
endif ()
if (EXISTS $ENV{TORCH_HOME})
    list(APPEND CMAKE_PREFIX_PATH $ENV{TORCH_HOME})
endif ()
if (CMAKE_PREFIX_PATH)
    list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
endif ()
message(STATUS "Looking for packages in : ${CMAKE_PREFIX_PATH}")
# ---

# check the compiler version
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (NOT CMAKE_CXX_COMPILER_VERSION) #Compiler version is not set on Ubuntu 12.02 (gcc 4.6)
        execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE CMAKE_CXX_COMPILER_VERSION)
    endif ()
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0)
        message("CXX COMPILER VERSION: ${CMAKE_CXX_COMPILER_ID} : ${CMAKE_CXX_COMPILER_VERSION}")
        message(FATAL_ERROR "Gadgetron requires GCC version >= 6.0")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.5)
        message(FATAL_ERROR "Gadgetron requires Clang version >= 3.5")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0)
        message(FATAL_ERROR "Gadgetron requires MSVC 2013 or above ")
    endif ()
else ()
    message(WARNING "Unsupported compiler!")
endif ()

#VERSIONING AND INSTALL PATHS
set(GADGETRON_VERSION_MAJOR 4)
set(GADGETRON_VERSION_MINOR 1)
set(GADGETRON_VERSION_PATCH 1)
set(GADGETRON_VERSION_STRING ${GADGETRON_VERSION_MAJOR}.${GADGETRON_VERSION_MINOR}.${GADGETRON_VERSION_PATCH})
set(GADGETRON_SOVERSION ${GADGETRON_VERSION_MAJOR}.${GADGETRON_VERSION_MINOR})
find_package(Git)
if (GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD WORKING_DIRECTORY
            ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE GADGETRON_GIT_SHA1 ERROR_VARIABLE GADGETRON_GIT_STDERR)
    string(STRIP "${GADGETRON_GIT_SHA1}" GADGETRON_GIT_SHA1)
    string(LENGTH "${GADGETRON_GIT_SHA1}" GADGETRON_GIT_SHA1_LEN)
    if (${GADGETRON_GIT_SHA1_LEN} LESS 40)
        message(WARNING "Could not determine SHA-1 hash: ${GADGETRON_GIT_STDERR}")
        set(GADGETRON_GIT_SHA1 "NA")
    endif ()
else ()
    set(GADGETRON_GIT_SHA1 "NA")
endif ()
set(GADGETRON_INSTALL_CMAKE_PATH share/gadgetron/cmake)
set(GADGETRON_INSTALL_CONFIG_PATH share/gadgetron/config)
set(GADGETRON_INSTALL_PYTHON_MODULE_PATH share/gadgetron/python)
set(GADGETRON_INSTALL_SCHEMA_PATH share/gadgetron/schema)
set(GADGETRON_INSTALL_INCLUDE_PATH include/gadgetron)
set(GADGETRON_INSTALL_CHROOT_SCRIPTS_PATH share/gadgetron/chroot)
if (WIN32)
    set(GADGETRON_VAR_DIR $ENV{APPDATA}/gadgetron CACHE PATH "Path where Gadgetron will store calibration data, etc" FORCE)
else ()
    set(GADGETRON_VAR_DIR /var/lib/gadgetron/ CACHE PATH "Path where Gadgetron will store calibration data, etc" FORCE)
endif ()
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# create the Gadgetron version header
configure_file("${CMAKE_SOURCE_DIR}/cmake/gadgetron_sha1.h.in" ${PROJECT_BINARY_DIR}/include/gadgetron_sha1.h @ONLY)
include_directories(${PROJECT_BINARY_DIR}/include)
install(FILES ${PROJECT_BINARY_DIR}/include/gadgetron_sha1.h DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

#Set the build type to Release if not specified
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
            FORCE)
endif ()

# build options for 64 bits system
if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    message("64bit system is found")
    set(HAS_64_BIT On CACHE BOOL "64bit build")
else ()
    message("32bit system is found")
    set(HAS_64_BIT Off CACHE BOOL "64bit build")
endif ()


#set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH ".:${CMAKE_INSTALL_PREFIX}/lib")

# whether to install dependencies
option(GADGETRON_INSTALL_DEPENDENCIES "Install gadgetron dependencies" Off)

# build options for OpenMP support
find_package(OpenMP)
option(USE_OPENMP "Use OpenMP" On)
if (OPENMP_FOUND)
    if (USE_OPENMP)
        message("OpenMP multithreading enabled")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        add_definitions(-DUSE_OMP)
    else ()
        message("OpenMP multithreading is supported, but disabled")
    endif ()
else ()
    message("OpenMP multithreading not supported")
endif ()

if (WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    add_definitions(-DWIN32 -D_WIN32 -D_WINDOWS -DWIN -D_AMD64_)
    add_definitions(-Zc:__cplusplus)
    #add_definitions(-DUNICODE -D_UNICODE)
    add_definitions("/wd4251") #disable warnings, 4251: needs to have dll-interface to be used by clients
    if (HAS_64_BIT)
        add_definitions(-DWIN64 -D_WIN64)
    endif ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /MP /bigobj /std:c++17")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W1")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /EHsc /MP /bigobj /w1")
    set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} /INCREMENTAL:NO")
    set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} /INCREMENTAL:NO")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
    set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
    set(CMAKE_STATIC_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
    set(CMAKE_MODULE_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
    # The two flags below is to fix Windows problems in relation to multiple defined operators new/delete and some constructors that are defined in our headers
    #  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /FORCE:MULTIPLE")
    #  set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /FORCE:MULTIPLE")
else ()
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    include(CheckCXXCompilerFlag)
    set(CMAKE_CXX_STANDARD 14)
    if (APPLE)
        if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        endif ()
    endif ()
    set(CXX_STANDARD 14)


    check_cxx_compiler_flag(-std=c++${CXX_STANDARD} COMPILER_SUPPORTS_stdcxx${CXX_STANDARD})
    if (COMPILER_SUPPORTS_stdcxx${CXX_STANDARD})
        set(CXX_STANDARD_FOUND 1)
        message(STATUS "Compiling using C++${CXX_STANDARD} standard")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${CXX_STANDARD}")
    endif ()
    check_cxx_compiler_flag(-Werror=return-type COMPILER_SUPPORTS_WRETURNTYPE)
    if (COMPILER_SUPPORTS_WRETURNTYPE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=return-type")
    endif ()
    check_cxx_compiler_flag(-Werror=terminate COMPILER_SUPPORTS_WTERMINATE)
    if (COMPILER_SUPPORTS_WTERMINATE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=terminate")
    endif ()
    check_cxx_compiler_flag(-Werror=non-virtual-dtor COMPILER_SUPPORT_WNON_VIRTUAL)
    if (COMPILER_SUPPORT_WNON_VIRTUAL)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=non-virtual-dtor")
    endif ()
endif ()

# whether to suppress compilation warnings
option(BUILD_SUPPRESS_WARNINGS "Build package while suppressing warnings" Off)
if (BUILD_SUPPRESS_WARNINGS)
    if (WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
    elseif (WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
    endif ()
endif ()

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
# necessary for Windows and RHEL <=6 systems
set(Boost_NO_BOOST_CMAKE ON)

if(WIN32)
add_definitions( -DHAVE_SNPRINTF) #Python decides to overwrite snprintf if this is not defined, because ofcourse it does.
endif()
# We actually only use system and thread explicitly, but they require linking in date_time and chrono
if (WIN32)
find_package(Boost 1.65.0 COMPONENTS system date_time chrono program_options filesystem timer exception REQUIRED )
else()
find_package(Boost 1.65.0 COMPONENTS system program_options filesystem timer REQUIRED )
endif()

add_library(boost INTERFACE)
set_target_properties(boost PROPERTIES
INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
message(${Boost_LIBRARIES})
#set_property(TARGET boost PROPERTY INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
target_link_libraries(boost INTERFACE ${Boost_LIBRARIES} ${CMAKE_DL_LIBS})
option(USE_CUDA "Build with CUDA support" On)

if (USE_CUDA)
find_package(CUDA_advanced)
endif()
if (CUDA_FOUND AND USE_CUDA)

if (NOT CUDA_cublas_device_LIBRARY)
list(REMOVE_ITEM CUDA_CUBLAS_LIBRARIES CUDA_cublas_device_LIBRARY-NOTFOUND)
endif()

add_definitions(-DUSE_CUDA)
set( GADGETRON_CUDA_FOUND_BOOL 1 )
include_directories( ${CUDA_INCLUDE_DIRS} )
#  if (NOT WIN32)
#    set(CUDA_PROPAGATE_HOST_FLAGS OFF)
#  endif ()
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} )
set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
#set(CUDA_VERBOSE_BUILD ON)

# Compile kernels for compute models 1.0 and 2.0 as default for Cuda 4.1
# Support compute model 3.0 from Cuda 4.2 and up
# Support compute model 3.5 from Cuda 5 and up

option(GADGETRON_CUDA_ALL_COMPUTE_MODEL "Build CUDA components for all computing models" OFF)
if (GADGETRON_CUDA_ALL_COMPUTE_MODEL)
message("Compiling CUDA components to support all compute model 3.5, 5.0, 5.2, 6.0, 6.1, 7.0")
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_35,code=sm_35")
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_50,code=sm_50")
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_52,code=sm_52")
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_60,code=sm_60")
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_61,code=sm_61")
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_70,code=sm_70")
endif ()

message("Compiling with ${CUDA_NVCC_FLAGS}")
else ()
message("CUDA not found. CUDA components will not be compiled.")
set( GADGETRON_CUDA_FOUND_BOOL 0 )
set( CUDA_NVCC_FLAGS "NA")
endif ()

find_package(GTest)

find_package(Armadillo REQUIRED )
add_library(armadillo INTERFACE)
set_property(TARGET armadillo PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${ARMADILLO_INCLUDE_DIR})
add_definitions(-DARMA_DONT_USE_WRAPPER)

find_package(PugiXML REQUIRED)
find_package(ISMRMRD REQUIRED)

link_directories(${ISMRMRD_LIBRARY_DIRS})
if (CUDA_FOUND AND USE_CUDA)
# bart needs cuda
    find_package(BART)
endif()

find_package(FFTW3 REQUIRED COMPONENTS single double)
add_library(FFTW INTERFACE)
target_link_libraries(FFTW INTERFACE ${FFTW3_LIBRARIES})
set_property(TARGET FFTW PROPERTY INTERFACE_LINK_DIRECTORIES ${FFTW3_INCLUDE_DIR})

set(GADGETRON_BUILD_RPATH "")

option(USE_MKL OFF)
if (USE_MKL)
    find_package(MKL REQUIRED)
    add_library(BLAS INTERFACE IMPORTED)
    set_property(TARGET BLAS PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${MKL_INCLUDE_DIR})
    set_property(TARGET BLAS PROPERTY INTERFACE_LINK_DIRECTORIES ${MKL_LIB_DIR})
    set_property(TARGET BLAS PROPERTY INTERFACE_LINK_LIBRARIES ${MKL_LIBRARIES})

    set(GADGETRON_BUILD_RPATH "${MKL_LIB_DIR};${MKL_COMPILER_LIB_DIR}")


  add_definitions(-DUSE_MKL)
  else()

  find_package(CBLAS REQUIRED)
	find_package(LAPACK)
	if (LAPACK_FOUND)
		message("LAPACK Found")
	else()
		find_package(clapack CONFIG REQUIRED)
		set(LAPACK_LIBRARIES f2c lapack)
	endif ()


endif()

if(EXISTS $ENV{TORCH_HOME})
    find_package(Torch)
    if (TORCH_FOUND)
        add_definitions(-DUSE_TORCH)
        message("Find Pytorch c++ library at ${TORCH_LIBRARIES}")
    else()
        message("Do not find Pytorch c++ library")
    endif ()
else()
    message("Do not find TORCH_HOME")
endif()

option(BUILD_PYTHON_SUPPORT "Build Python gadgets etc" On)

if (BUILD_PYTHON_SUPPORT)

find_package(PythonLibs 3  REQUIRED)

find_package(NumPy REQUIRED)
message("Found numpy : ${NUMPY_INCLUDE_DIRS} ${NUMPY_LIBRARIES} ${PYTHON_INCLUDE_DIRS} ${PYTHON_LIBRARIES}")


add_library(python INTERFACE)
set_property(TARGET python PROPERTY
INTERFACE_INCLUDE_DIRECTORIES ${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIRS})
message("Boost_VERSION ${Boost_VERSION_STRING}")
if (Boost_VERSION_STRING VERSION_LESS 1.67.0)
    find_package(Boost 1.65.0 COMPONENTS python3 REQUIRED)
    else()
    string(REGEX MATCH "^3\\.([0-9]+)\\.[0-9]+" PYTHON_MINOR_VERSION ${PYTHONLIBS_VERSION_STRING} )
    set(PYTHON_MINOR_VERSION ${CMAKE_MATCH_1})
    find_package(Boost 1.65.0 COMPONENTS "python3${PYTHON_MINOR_VERSION}" REQUIRED)
endif()
set_property(TARGET python PROPERTY
INTERFACE_INCLUDE_DIRECTORIES ${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIRS})
target_link_libraries(python INTERFACE ${PYTHON_LIBRARIES} ${NUMPY_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(boost INTERFACE python)

endif()

message("Searching for OpenGL, GLEW, GLUT, and Qt. These libraries are only used in a single standalone application and are thus non-essential.")
if(WIN32)
message("For Windows users in particular, for ease of installation we do not reccomend installing these libraries.")
endif()

find_package(OpenGL)
find_package(GLEW)
find_package(GLUT)
find_package(Qt4 4.6)
find_package(PLplot)

add_subdirectory(core)
add_subdirectory(toolboxes)
add_subdirectory(gadgets)
add_subdirectory(apps)

option(BUILD_WITH_GPERFTOOLS_PROFILER "Build with the gperftools profiler." Off)
if (BUILD_WITH_GPERFTOOLS_PROFILER)
message("Adding gpertools to build.")
find_package(Gperftools REQUIRED)
endif()

add_subdirectory(test)
add_subdirectory(cmake)
add_subdirectory(chroot)

option(BUILD_DOCUMENTATION "Build Gadgetron documentation." Off)
if(BUILD_DOCUMENTATION)
    add_subdirectory(doc)
endif()

# install dependencies for WIN32
if (WIN32)
if (GADGETRON_INSTALL_DEPENDENCIES)
include(${CMAKE_SOURCE_DIR}/cmake/InstallWinDependencies.cmake)
endif ()
endif ()

if (UNIX)
if (NOT APPLE)
if (GADGETRON_INSTALL_DEPENDENCIES)
include(${CMAKE_SOURCE_DIR}/cmake/InstallLinuxDependencies.cmake)
endif ()
endif ()
endif ()

#  ---   Main Library  (end) ----

# Create package
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)
include(${CMAKE_SOURCE_DIR}/cmake/gadgetron_cpack.cmake)
if(CPACK_GENERATOR)
message(STATUS "Found CPack generators: ${CPACK_GENERATOR}")
configure_file("${CMAKE_SOURCE_DIR}/cmake/cpack_options.cmake.in" ${GADGETRON_CPACK_CFG_FILE} @ONLY)
set(CPACK_PROJECT_CONFIG_FILE ${GADGETRON_CPACK_CFG_FILE})
include(CPack)
endif()
