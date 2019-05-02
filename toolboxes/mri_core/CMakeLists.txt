if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_MRI_CORE__)
endif ()

include_directories(
    ${Boost_INCLUDE_DIR} 
    ${ARMADILLO_INCLUDE_DIRS} 
    ${ACE_INCLUDE_DIR} 
    ${ISMRMRD_INCLUDE_DIR}
    ${FFTW3_INCLUDE_DIR}
    ${MKL_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/apps/gadgetron
    ${CMAKE_BINARY_DIR}/apps/gadgetron
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_BINARY_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/klt/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/operators
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math 
    ${CMAKE_SOURCE_DIR}/toolboxes/image_io
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
)

set( mri_core_header_files
        mri_core_export.h
        mri_core_def.h
        mri_core_data.h
        mri_core_utility.h
        mri_core_kspace_filter.h
        mri_core_grappa.h 
        mri_core_spirit.h 
        mri_core_coil_map_estimation.h 
        mri_core_dependencies.h 
        mri_core_acquisition_bucket.h
        mri_core_girf_correction.h
        mri_core_partial_fourier.h )

set( mri_core_source_files
        mri_core_utility.cpp 
        mri_core_grappa.cpp 
        mri_core_spirit.cpp 
        mri_core_kspace_filter.cpp
        mri_core_coil_map_estimation.cpp 
        mri_core_dependencies.cpp 
        mri_core_partial_fourier.cpp mri_core_girf_correction.cpp )

add_library(gadgetron_toolbox_mri_core SHARED 
     ${mri_core_header_files} ${mri_core_source_files} )

set_target_properties(gadgetron_toolbox_mri_core PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
set_target_properties(gadgetron_toolbox_mri_core PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(gadgetron_toolbox_mri_core 
                    gadgetron_toolbox_cpucore 
                    gadgetron_toolbox_cpucore_math
                    ${ARMADILLO_LIBRARIES} 
                    gadgetron_toolbox_cpufft 
                    gadgetron_toolbox_cpuklt 
                    gadgetron_toolbox_image_analyze_io 
                    )

install(TARGETS gadgetron_toolbox_mri_core DESTINATION lib COMPONENT main)

install(FILES
  ${mri_core_header_files}
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

# if there is python, add python binding
if (Boost_PYTHON3_FOUND AND PYTHONLIBS_FOUND AND NUMPY_FOUND)

    include_directories(
      ${CMAKE_SOURCE_DIR}/toolboxes/python
      ${PYTHON_INCLUDE_PATH}
      ${NUMPY_INCLUDE_DIRS}
      ${Boost_INCLUDE_DIR}
      )

      add_library(gadgetron_toolbox_mri_core_python MODULE 
                mri_core_grappa_python.h
                mri_core_grappa_python.cpp
                mri_core_python.cpp )

      target_link_libraries(gadgetron_toolbox_mri_core_python
                gadgetron_toolbox_cpucore 
                gadgetron_toolbox_cpucore_math
                ${ARMADILLO_LIBRARIES} 
                gadgetron_toolbox_cpufft 
                gadgetron_toolbox_cpuklt 
                gadgetron_toolbox_image_analyze_io 
                gadgetron_toolbox_python 
                gadgetron_toolbox_mri_core 
                ${PYTHON_LIBRARIES}
                ${Boost_LIBRARIES}
                ${MKL_LIBRARIES}
                )

    if (WIN32)
        set_target_properties(gadgetron_toolbox_mri_core_python PROPERTIES SUFFIX .pyd)
    endif()

    set_target_properties(gadgetron_toolbox_mri_core_python PROPERTIES PREFIX "")

    install(TARGETS gadgetron_toolbox_mri_core_python DESTINATION ${GADGETRON_INSTALL_PYTHON_MODULE_PATH} COMPONENT main)

endif ()

