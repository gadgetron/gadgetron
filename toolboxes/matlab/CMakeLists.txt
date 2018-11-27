
include_directories(${MATLAB_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/toolboxes/mri_core 
${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math 
${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image 
${CMAKE_SOURCE_DIR}/toolboxes/core/
${CMAKE_BINARY_DIR}/toolboxes/core/
${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
${CMAKE_SOURCE_DIR}/toolboxes/core/gpu 
${CMAKE_SOURCE_DIR}/toolboxes/operators
${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu 
${ISMRMRD_INCLUDE_DIR}
${Boost_INCLUDE_DIR}
${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/ )



add_library(gadgetron_toolbox_matlab SHARED  matlab_export.h MatlabUtils.h MatlabUtils.cpp )


set_target_properties(gadgetron_toolbox_matlab PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})



target_link_libraries(
    gadgetron_toolbox_matlab
    gadgetron_toolbox_log
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_cpucore_math
    ${MATLAB_LIBRARIES}
    ${ISMRMRD_LIBRARIES}

)

if (UNIX AND CUDA_FOUND)
#add_subdirectory(@gadgetron)
endif()

install(TARGETS gadgetron_toolbox_matlab DESTINATION lib COMPONENT main)

install(FILES MatlabUtils.h DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
